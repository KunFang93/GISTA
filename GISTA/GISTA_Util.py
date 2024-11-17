import sys
import time
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from collections import Counter, defaultdict
from scipy import stats
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from scipy.stats import pearsonr
import logomaker
from functools import reduce

# v5 add two conditions comparison
def generate_comparison_dict(samples_df, comparison_string):
    # Group samples by GroupID
    group_samples = samples_df.groupby("GroupID")["IndvidualID"].apply(list).to_dict()
    # Parse comparison string and create the dictionary
    comparison_dict = {}
    comparison_types = {}
    combined_samples = {}
    # Add intra-group comparisons for all groups
    for group, individuals in group_samples.items():
        comparison_dict[f"Intra-{group}"] = [individuals, individuals]
    for pair in comparison_string.split(";"):
        group, types = pair.split(",")
        cur_comparison = []
        # Handle combined groups in "group" and "types"
        if "-" in group:
            combined_name = group.replace("-", "")
            combined_samples[combined_name] = [
                individual
                for subtype in group.split("-")
                for individual in group_samples[subtype.strip()]
            ]
            comparison_dict[f"Intra-{combined_name}"] = [
                combined_samples[combined_name],
                combined_samples[combined_name],
            ]
            if f"Intra-{combined_name}" not in cur_comparison:
                cur_comparison.append(f"Intra-{combined_name}")

        if "-" in types:
            combined_name = types.replace("-", "")
            combined_samples[combined_name] = [
                individual
                for subtype in types.split("-")
                for individual in group_samples[subtype.strip()]
            ]
            comparison_dict[f"Intra-{combined_name}"] = [
                combined_samples[combined_name],
                combined_samples[combined_name],
            ]
            if f"Intra-{combined_name}" not in cur_comparison:
                cur_comparison.append(f"Intra-{combined_name}")
        # Handle inter-group comparisons
        if "-" in group:
            group = group.replace("-", "")
            first_group = combined_samples[group.strip()]
        else:
            first_group = group_samples[group.strip()]
        if "-" in types:
            types = types.replace("-", "")
            second_group = combined_samples[types.strip()]
        else:
            second_group = group_samples[types.strip()]
        comparison_name = f"{group}vs{types}"
        comparison_dict[comparison_name] = [first_group, second_group]
        cur_comparison.insert(0, comparison_name)  # Insert inter-group comparison at the start
        # Add intra-group comparisons for individual groups in "group" and "types"
        for grp in group.split("-") + types.split("-"):
            grp = grp.strip()
            intra_name = f"Intra-{grp}"
            if intra_name not in cur_comparison:
                cur_comparison.append(intra_name)
        comparison_types[comparison_name] = cur_comparison
    return comparison_dict, comparison_types


class TadsArray(object):
    def __init__(self, tads_files, samples, binsize):
        # it is important that tads_files and samples in the same order
        # flexibility represents how many bins that allow to consider two ids are closed. e.g. if flexibilty = 1, then 17 and 18
        # are considered as closed but not for 16 and 18
        self.tads_files = tads_files
        self.samples = samples
        self.binsize = binsize
    def loading_files(self):
        print("Loading files")
        data_tads = {sample:0 for sample in self.samples}
        tad_ids = []
        for idx,file in enumerate(self.tads_files):
            current_df = pd.read_csv(file,sep='\t')
            # drop row that has same from.coord and to.coord, annoying topdom bugs
            current_df = current_df.drop(current_df.loc[current_df.loc[:,'from.coord']==current_df.loc[:,'to.coord'],:].index.values)
            current_df.reset_index(inplace=True)
            # make sure same data type of "to.id"
            current_df["to.id"] = pd.to_numeric(current_df["to.id"], downcast="float")
            combined_mark = list(current_df.loc[:,'chr']+'_'+current_df.loc[:,'to.id'].astype(str))
            tad_ids.append(combined_mark)
            current_df.loc[:,'mark'] = combined_mark
            # assign the vector value for each region
            current_df.loc[:,'vec'] = current_df.loc[:,'size']/self.binsize
            current_df.loc[current_df.loc[:,'tag']=='gap','vec'] = -current_df.loc[current_df.loc[:,'tag']=='gap','vec']
            current_df.loc[current_df.loc[:,'tag']=='boundary','vec'] = -current_df.loc[current_df.loc[:,'tag']=='boundary','vec']
            data_tads[self.samples[idx]] = current_df
        return data_tads, tad_ids
    def tads_basic_stats(self, data_tads):
        data_stat_tads = {'Sample':[],'domain':[],'boundary':[],'gap':[],'mean.tads.size':[]}
        data_dist_tads = {'Sample':[],'size':[]}
        for sample in self.samples:
            current_df = data_tads[sample]
            # get some statistics
            current_size = current_df.loc[current_df.loc[:,'tag']=='domain','size']
            data_dist_tads['Sample'] += [sample] * len(current_size)
            data_dist_tads['size'] += list(current_size)
            data_stat_tads['Sample'].append(sample)
            data_stat_tads['mean.tads.size'].append(current_size.mean())
            for (k,v) in Counter(current_df['tag']).items():
                data_stat_tads[k].append(v)
        data_dist_tads_df = pd.DataFrame.from_dict(data_dist_tads)
        data_stat_tads_df = pd.DataFrame.from_dict(data_stat_tads)
        return data_dist_tads_df,data_stat_tads_df
    def _split_by_common(self, data_tads, common_coor):
        # return data_tads_split (data_tads_split['NT1']['chr1'][0])
        total_num_tads_array = {sample:0 for sample in self.samples}
        data_tads_split ={sample:{} for sample in self.samples}
        for sample in self.samples:
            current_df = data_tads[sample]
            # group by chromosomes
            current_gb = current_df.groupby('chr')
            for chrom in current_gb.groups:
                current_chr_df = current_df.loc[current_gb.groups[chrom],:]
                # +1 is essential for the correct splits
                idx_list = current_chr_df[current_chr_df.loc[:,'mark'].isin(common_coor)].index.values + 1
                idx_list_sort = sorted(idx_list)
                idx_list_sort = idx_list_sort - current_chr_df.index.values[0]
                # avoid generating empty array
                if idx_list_sort[-1] == len(current_chr_df):
                    idx_list_sort = idx_list_sort[:-1]
                # merge second last tads array and last tads array for each chromosome (empirical)
                idx_list_sort = idx_list_sort[:-1]
                data_tads_split[sample][chrom] = np.split(current_chr_df,idx_list_sort)
                total_num_tads_array[sample] += len(data_tads_split[sample][chrom])
                # make sure each chromosome in data_tads_split[sample] has same number of tads array
        test_total = list(total_num_tads_array.values())[0]
        for sample in self.samples:
            if total_num_tads_array[sample] != test_total:
                print("Encounter different number of tads array, exit")
                exit(1)
        return data_tads_split
    def _tads_matrix(self,data_tads_split,chroms):
        vector_matrices = {}
        # vector_annot['tads_id'] = [chr,start,end,{'NT1':"d12,b10,g1,d20"}]
        vector_annot = defaultdict(list)
        tads_id_list = []
        for chrom in chroms:
            for i in range(len(data_tads_split[self.samples[0]][chrom])):
                vec_list = []
                try:
                    tads_id = chrom+'_'+str(i)
                    tads_id_list.append(tads_id)
                    # chr
                    vector_annot[tads_id].append(data_tads_split[self.samples[0]][chrom][i].iloc[0,1])
                    # start
                    vector_annot[tads_id].append(data_tads_split[self.samples[0]][chrom][i].iloc[0,3])
                    # end
                    vector_annot[tads_id].append(data_tads_split[self.samples[0]][chrom][i].iloc[-1,5])
                    current_sample_annot = {}
                    for sample in self.samples:
                        current_annot_list = list(data_tads_split[sample][chrom][i].loc[:,'tag'].str[0] +
                                                  (data_tads_split[sample][chrom][i].loc[:,'size']/binsize).astype(str))
                        current_sample_annot[sample]= ','.join(current_annot_list)
                        current_vec = data_tads_split[sample][chrom][i].loc[:,'vec']
                        current_vec.index.values[:] = np.arange(len(current_vec))
                        vec_list.append(data_tads_split[sample][chrom][i].loc[:,'vec'])
                    vector_annot[tads_id].append(current_sample_annot)
                    current_vec_matrix = pd.concat(vec_list,axis=1)
                    current_vec_matrix.columns = self.samples
                    vector_matrices[chrom+'_'+str(i)] = current_vec_matrix.T
                except IndexError:
                    # last element is blank
                    print("Skip blank elements")
                    continue
        return vector_matrices, vector_annot, tads_id_list
    def _matrix2array(self,matrix):
        matrix_fill = matrix.fillna(0)
        matrix_dict = matrix_fill.T.to_dict('list')
        final_dict = {}
        for key in matrix_dict:
            final_dict[key] = [i for i in matrix_dict[key] if i !=0]
        return final_dict
    def _softseg(self,trt_list,ctrl_list):
        # flexibility is the parameter define the degree of the softness of the border
        flexibility = 2
        trt_abs = np.abs(trt_list)
        ctrl_abs = np.abs(ctrl_list)
        trt_cumsum = np.cumsum(trt_abs)
        ctrl_cumsum = np.cumsum(ctrl_abs)
        pile_trt = np.tile(trt_cumsum,(len(ctrl_cumsum),1))
        index_trt = pile_trt - np.asarray(ctrl_cumsum)[:,None]
        index_trt_abs = np.abs(index_trt)
        # find the index of the soft border
        ctrl_idx,trt_idx = np.where(index_trt_abs <= flexibility)
        # print(trt_idx,ctrl_idx)
        # avoid consecutive soft border
        new_ctrl_idx = []
        new_trt_idx = []
        for idx,value in enumerate(ctrl_idx):
            if value not in new_ctrl_idx and trt_idx[idx] not in new_trt_idx:
                new_ctrl_idx.append(value)
                new_trt_idx.append(trt_idx[idx])
            else:
                continue
        # print(new_trt_idx,new_ctrl_idx)
        new_ctrl_idx = np.array(new_ctrl_idx)
        new_trt_idx = np.array(new_trt_idx)
        new_ctrl_idx += 1
        new_trt_idx += 1
        new_ctrl_idx = np.concatenate([[0],new_ctrl_idx])
        new_trt_idx = np.concatenate([[0],new_trt_idx])
        ctrl_tads_num = np.diff(new_ctrl_idx)
        trt_tads_num = np.diff(new_trt_idx)
        ctrl_tads_split_idx = np.cumsum(ctrl_tads_num)
        trt_tads_split_idx = np.cumsum(trt_tads_num)
        return trt_tads_split_idx[:-1],ctrl_tads_split_idx[:-1]
    def new_idx(self,split_idx_dict,minlen_sample):
        final_split_idx = {sample:0 for sample in self.samples}
        seg_idx = [split_idx_dict[minlen_sample][sample][0] for sample in split_idx_dict[minlen_sample]]
        # find the common seg for all pairwised comparison
        common_seg_idx = sorted(list(set(seg_idx[0]).intersection(*seg_idx)))
        final_split_idx[minlen_sample] = common_seg_idx
        for sample in split_idx_dict[minlen_sample]:
            target_seg_idx = split_idx_dict[minlen_sample][sample][1]
            current_seg_idx =  split_idx_dict[minlen_sample][sample][0]
            new_idx = [current_seg_idx.tolist().index(idx) for idx in common_seg_idx]
            final_split_idx[sample] = [target_seg_idx[idx] for idx in new_idx]
        return final_split_idx
    def find_subcommon(self,input_dict,samples):
        split_idx_dict = {sample:{} for sample in self.samples}
        min_len = len(min(input_dict.items(),key=lambda x:len(x[1]))[1])
        minlen_samples = [sample for sample in samples if len(input_dict[sample])==min_len]
        final_seg = {}
        # store all pairwise segmentation index
        for trt_ctrl_sample in itertools.combinations(samples,2):
            trt_sample = trt_ctrl_sample[0]
            ctrl_sample = trt_ctrl_sample[1]
            trt_split_idx, ctrl_split_idx = self._softseg(input_dict[trt_sample],input_dict[ctrl_sample])
            if len(trt_split_idx)==0 or len(ctrl_split_idx) ==0:
                return {sample:[np.array(input_dict[sample])] for sample in input_dict}
            else:
                split_idx_dict[trt_sample][ctrl_sample] = [trt_split_idx,ctrl_split_idx]
                split_idx_dict[ctrl_sample][trt_sample] = [ctrl_split_idx,trt_split_idx]
        if len(minlen_samples) ==1:
            # find the shortest segment
            final_split_idx = self.new_idx(split_idx_dict,minlen_samples[0])
        else:
            idx_len = []
            tmp_store = {}
            for idx,minlen_sample in enumerate(minlen_samples):
                split_idx = self.new_idx(split_idx_dict,minlen_sample)
                idx_len.append(len(split_idx[minlen_sample]))
                tmp_store[idx] = split_idx
            minidx = int(np.argmin(idx_len))
            final_split_idx = tmp_store[minidx]
        for sample in samples:
            final_seg[sample] = np.split(input_dict[sample],final_split_idx[sample])
        return final_seg
    def tads_array_chunk(self, data_tads, tad_ids):
        # find common coordinates id for all samples
        common_coor = list(set(tad_ids[0]).intersection(*tad_ids))
        print("The number of common coordinates(0 flexiblity): {}".format(len(common_coor)))
        chroms = [i for i in data_tads[self.samples[0]].groupby('chr').groups]
        print("First round segmentation")
        data_tads_split = self._split_by_common(data_tads,common_coor)
        # construct tad vector matrix, each row represent the tads array for a sample (length might be varied)
        print("Construct vector matrices")
        vector_matrices, vector_annot, tads_id_list = self._tads_matrix(data_tads_split,chroms)
        print("Second round segmentation")
        tads_sub_array_dict = {}
        tads_sub_annot_dict = {}
        tads_sub_id_list =[]
        for array_mat_id in tqdm(vector_matrices):
            current_df = vector_matrices[array_mat_id]
            current_dict = self._matrix2array(current_df)
            current_split_seg = self.find_subcommon(current_dict,samples)
            chunk_num = len(current_split_seg[self.samples[0]])
            for idx in range(chunk_num):
                tads_sub_id = "{}_{}".format(array_mat_id,idx)
                tads_sub_id_list.append(tads_sub_id)
                tads_sub_array_dict[tads_sub_id] = {sample:current_split_seg[sample][idx] for sample in self.samples}
                if idx != 0:
                    tads_sub_annot_dict[tads_sub_id] = {sample:[vector_annot[array_mat_id][0],
                                                                vector_annot[array_mat_id][1] + self.binsize * np.sum(np.abs(np.concatenate(current_split_seg[sample][:idx]))),
                                                                vector_annot[array_mat_id][1] + self.binsize * np.sum(np.abs(np.concatenate(current_split_seg[sample][:idx+1])))]
                                                        for sample in self.samples}
                else:
                    tads_sub_annot_dict[tads_sub_id] = {sample:[vector_annot[array_mat_id][0],
                                                                vector_annot[array_mat_id][1],
                                                                vector_annot[array_mat_id][1] + self.binsize * np.sum(np.abs(current_split_seg[sample][0]))]
                                                        for sample in self.samples}
        return tads_sub_array_dict, tads_sub_annot_dict, tads_sub_id_list

class TadsMatrix(object):
    def __init__(self,tads_array_dict,samples,tads_id_list,onerep=False, pseudorep =5):
        self.flexibility = 2
        self.tads_id_list = tads_id_list
        if not onerep:
            self.tads_array_dict = tads_array_dict
            self.samples = samples
        else:
            # make identical pseudo-rep
            new_tads_array_dict = {tads_id:{} for tads_id in tads_id_list}
            new_samples = [f'{sample}{i}' for sample in samples for i in range(pseudorep)]
            for tads_id in tads_id_list:
                cur_dict = tads_array_dict[tads_id]
                for key in cur_dict:
                    for i in range(pseudorep):
                        new_tads_array_dict[tads_id][f'{key}{i}'] = cur_dict[key]
            self.tads_array_dict = new_tads_array_dict
            self.samples = new_samples
    def _neo_del_score(self,input_dict):
        # neo positive; del negative
        sum_series = pd.DataFrame.from_dict(input_dict,orient='index').sum(axis=1)
        # row ctrl, column trt, mat_row_expand is dup sum_series in rows
        mat_row_expand = np.tile(sum_series.tolist(),(len(self.samples),1))
        mat_col_expand = mat_row_expand.T
        # find diff between any two samples
        nd_mat = pd.DataFrame((mat_row_expand - mat_col_expand)/2,index=self.samples,columns=self.samples)
        # if neo_del less than flexibility, make it 0. Might change to some relative conserve information in the future
        nd_mat[nd_mat.abs()<=self.flexibility] = 0
        return nd_mat
    def _split_fuse_score(self,trt_list,ctrl_list):
        if len(trt_list) > len(ctrl_list):
            split_fuse_score = sum(np.sort(np.abs(trt_list))[::-1][len(ctrl_list):])
        elif len(trt_list) < len(ctrl_list):
            split_fuse_score = -sum(np.sort(np.abs(ctrl_list))[::-1][len(trt_list):])
        else:
            split_fuse_score = 0.0
        if np.abs(split_fuse_score) <= self.flexibility:
            split_fuse_score = 0.0
        return split_fuse_score
    def matrix_build(self,input_dict):
        # return 'skew-symmetric' matrix with shape sample x sample. The 'rows' used as the ctrl, for exmaple: the element in ['NT1','NT2'] means
        # use 'NT1' as ctrl, the tads array change in 'NT2'; each element in the matrix is in two dimensional [neo/del,split/fuse].
        sf_df_raw = pd.DataFrame(0,index=self.samples,columns=self.samples)
        nd_df_raw = self._neo_del_score(input_dict)
        nd_df = nd_df_raw.abs()
        for sample in itertools.combinations(self.samples,2):
            trt_sample = sample[0]
            ctrl_sample = sample[1]
            trt_array = input_dict[trt_sample]
            ctrl_array = input_dict[ctrl_sample]
            split_fuse_score = self._split_fuse_score(trt_array,ctrl_array)
            sf_df_raw.loc[trt_sample,ctrl_sample] = -split_fuse_score
        sf_df_raw = pd.DataFrame(np.triu(sf_df_raw) - np.triu(sf_df_raw,1).T,index=self.samples,columns=self.samples)
        sf_df = sf_df_raw.abs()
        detail_info = {'nd_raw':nd_df_raw,'sf_raw':sf_df_raw}
        return nd_df, sf_df, detail_info
    def fit(self,comparison_dict):
        # out_matrix_dict = {}
        # stats_dict is feature dict
        stats_dict = {'tads_id':self.tads_id_list}
        for compare_type in comparison_dict:
            stats_type_nd = compare_type + '_nd'
            stats_type_sf = compare_type + '_sf'
            stats_dict[stats_type_nd] = []
            stats_dict[stats_type_sf] = []
        print("Building featured matrix")
        total_info = {'total_detail':{},'total_nd':{},'total_sf':{}}
        for tads_array_id in tqdm(self.tads_id_list):
            current_dict = self.tads_array_dict[tads_array_id]
            current_nd, current_sf, detail_info = self.matrix_build(current_dict)
            total_info['total_detail'][tads_array_id] = detail_info
            total_info['total_nd'][tads_array_id] = current_nd
            total_info['total_sf'][tads_array_id] = current_sf
            # fill the stats_dict
            for compare_type in comparison_dict:
                stats_type_nd = compare_type + '_nd'
                stats_type_sf = compare_type + '_sf'
                trt_samples = comparison_dict[compare_type][0]
                ctrl_samples = comparison_dict[compare_type][1]
                # inter and intra tissue type comparison
                current_part_nd = current_nd.loc[ctrl_samples,trt_samples]
                current_part_sf = current_sf.loc[ctrl_samples,trt_samples]
                if trt_samples != ctrl_samples:
                    # inter
                    stats_dict[stats_type_nd].append(np.round(current_part_nd.stack().mean(),3))
                    stats_dict[stats_type_sf].append(np.round(current_part_sf.stack().mean(),3))
                else:
                    tri_nd = np.array(current_part_nd)[np.triu_indices(len(trt_samples),k=1)]
                    tri_sf = np.array(current_part_sf)[np.triu_indices(len(trt_samples),k=1)]
                    stats_dict[stats_type_nd].append(np.round(np.mean(tri_nd),3))
                    stats_dict[stats_type_sf].append(np.round(np.mean(tri_sf),3))
        stats_df = pd.DataFrame(stats_dict)
        return total_info, stats_df

def permutation_test(df, gcomp, lvl1_coef=0.7, lvl2_coef=0.3, n_permutations=1000):
    """Performs a Monte Carlo permutation test on a row of data and returns the p-value."""
    data = df.copy()
    # Set the seed for the random number generator
    np.random.seed(42)
    lvl1 = gcomp.split('.')[2]
    lvl2 = gcomp.split('.')[0]
    # Q99 of TT.vs.NT_nd
    Qnd_group = data['{}_nd'.format(gcomp)].quantile(0.99)
    Qsf_group = data['{}_sf'.format(gcomp)].quantile(0.99)
    # Q99 value of Intra-NT_nd
    Qnd_lvl1 = data['Intra-{}_nd'.format(lvl1)].quantile(0.99)
    Qsf_lvl1 = data['Intra-{}_sf'.format(lvl1)].quantile(0.99)
    Qnd_lvl2 = data['Intra-{}_nd'.format(lvl2)].quantile(0.99)
    Qsf_lvl2 = data['Intra-{}_sf'.format(lvl2)].quantile(0.99)
    # print('Qnd_group:{},Qnd_lvl1:{},Qnd_lvl2:{}'.format(Qnd_group,Qnd_lvl1,Qnd_lvl2))
    # print('Qsf_group:{},Qsf_lvl1:{},Qsf_lvl2:{}'.format(Qsf_group, Qsf_lvl1, Qsf_lvl2))
    # Calculate the observed test statistic
    observed_statistic = (data['{}_nd'.format(gcomp)]/Qnd_group + data['{}_sf'.format(gcomp)]/Qsf_group) - \
                         lvl1_coef * (data['Intra-{}_nd'.format(lvl1)]/Qnd_lvl1 + data['Intra-{}_sf'.format(lvl1)]/Qsf_lvl1) - \
                         lvl2_coef * (data['Intra-{}_nd'.format(lvl2)]/Qnd_lvl2 + data['Intra-{}_sf'.format(lvl2)]/Qsf_lvl2)

    # Generate a distribution of test statistics under the null hypothesis
    permuted_obs = [pd.Series(np.random.permutation(observed_statistic)) for i in range(n_permutations)]
    permuted_obs_df = pd.concat(permuted_obs,axis=1)

    # Calculate the p-value as the proportion of permuted test statistics that are as or more extreme than the observed test statistic
    p_values = (permuted_obs_df >= observed_statistic.values[:, None]).mean(axis=1)
    return p_values

def cal_pval(feature_df, gcomp):
    max_d_number = 0
    final_pval = 0
    # optimize lvl1 and lvl2 coefficient to get the maximum number of dchange
    only_inter_pval = 0
    for i in range(10):
        for j in range(10):
            lvl1_coef = i / 10
            lvl2_coef = j / 10
            p_values = permutation_test(df=feature_df, lvl1_coef=lvl1_coef, lvl2_coef=lvl2_coef,
                                                    gcomp=gcomp)
            sig_num = len(p_values[p_values < 0.05])
            if i == 0 and j == 0:
                only_inter_pval = permutation_test(df=feature_df, lvl1_coef=lvl1_coef, lvl2_coef=lvl2_coef,
                                                    gcomp=gcomp)
            if max_d_number < sig_num:
                max_d_number = sig_num
                final_pval = p_values
            else:
                continue
    return [only_inter_pval,final_pval]

def class_process(input_df,group_level,individual_level1,individual_level2,group_diff_cut,individual_diff_cut):
    # group_level TT.vs.NT_nd
    df = input_df.copy()
    group_cut_up = input_df.stack()[input_df.stack()!=0].quantile(group_diff_cut)
    indi_cut_up = input_df.stack()[input_df.stack()!=0].quantile(individual_diff_cut)
    group_cut_down = input_df.stack()[input_df.stack()!=0].quantile(0.1)
    indi_cut_down = group_cut_down
    print("up_cut_val: {}; down_up_val: {}".format(group_cut_up,group_cut_down))
    df.loc[df.loc[:,group_level]>group_cut_down,group_level+'_class'] = 'medium'
    df.loc[df.loc[:,individual_level1]>indi_cut_down,individual_level1+'_class'] = 'medium'
    df.loc[df.loc[:,individual_level2]>indi_cut_down,individual_level2+'_class'] = 'medium'
    df.loc[df.loc[:,group_level]>group_cut_up,group_level+'_class'] = 'high'
    df.loc[df.loc[:,individual_level1]>indi_cut_up,individual_level1+'_class'] = 'high'
    df.loc[df.loc[:,individual_level2]>indi_cut_up,individual_level2+'_class'] = 'high'
    df = df.fillna('low')
    return df, [group_cut_up,group_cut_down]

def combine_nd_sf_annot(tads_id,nd_marked_df,sf_marked_df,group_level,individual_level1,individual_level2):
    # group_level TT.vs.NT
    # diff_cut 0.7 means we define top 30% as the drastically change
    combine_dict = {'high_high':'high','high_medium':'high','high_low':'high','medium_high':'high','medium_medium':'medium','medium_low':'medium',
                    'low_high':'high','low_medium':'medium','low_low':'low'}
    types_dict = {'high_high':'Mixed','high_medium':'Mixed','high_low':'ND','medium_high':'Mixed','medium_medium':'NotAvail','medium_low':'NotAvail',
                  'low_high':'SF','low_medium':'NotAvail','low_low':'NotAvail'}
    final_df = pd.DataFrame(tads_id)
    final_df[group_level] = (nd_marked_df['{}_nd_class'.format(group_level)] + '_' +
                             sf_marked_df['{}_sf_class'.format(group_level)]).map(combine_dict)
    final_df[individual_level1] = (nd_marked_df['{}_nd_class'.format(individual_level1)] + '_' +
                                   sf_marked_df['{}_sf_class'.format(individual_level1)]).map(combine_dict)
    final_df[individual_level2] = (nd_marked_df['{}_nd_class'.format(individual_level2)] + '_' +
                                   sf_marked_df['{}_sf_class'.format(individual_level2)]).map(combine_dict)
    final_df['GroupChange'] = 'MV'
    final_df.loc[final_df[group_level]=='high','GroupChange'] = 'SV'
    final_df.loc[(final_df[group_level]=='low') & (final_df[individual_level1]!='high') & (final_df[individual_level2]!='high'),'GroupChange'] = 'C'
    # add ND, SF and ND_SF_Mixed
    final_df['ChangeTypes'] = (nd_marked_df['{}_nd_class'.format(group_level)] + '_' +
                               sf_marked_df['{}_sf_class'.format(group_level)]).map(types_dict)
    # add Individual-high, Individual-M & C
    final_df['IndividualChange'] = 'NotAvail'
    final_df.loc[(final_df[group_level] == 'high') & (final_df[individual_level2] == 'high'),'IndividualChange'] = 'HS'
    final_df.loc[
        (final_df[group_level] == 'high') & (final_df[individual_level2] != 'high'), 'IndividualChange'] = 'LS'
    # final_df.loc[(final_df[group_level]=='low') & (final_df[individual_level1]='low') & (final_df[individual_level2]=='low'),'bio_mark'] = 'C'
    # final_df.loc[(final_df[group_level]=='low') & (final_df[individual_level1]=='meidum') & (final_df[individual_level2]=='low'),'bio_mark'] = 'C'
    return final_df

def coor_annot(features_df_marked, tads_sub_id_list,tads_sub_annot_dict):
    cols = list(features_df_marked.columns)
    tads_subarray_start = []
    tads_subarray_end = []
    tads_subarray_chr = []
    for tads_sub_id in tads_sub_id_list:
        # to avoid overlapping, pick the largest start and smallest end coodinates
        min_start = min(tads_sub_annot_dict[tads_sub_id].items(),key=lambda x:x[1][1])[1][1]
        max_end = max(tads_sub_annot_dict[tads_sub_id].items(),key=lambda x:x[1][2])[1][2]
        tads_subarray_start.append(int(min_start))
        tads_subarray_end.append(int(max_end))
        tads_subarray_chr.append(tads_sub_annot_dict[tads_sub_id][samples[0]][0])
    features_df_marked.loc[:,'start(min)'] = tads_subarray_start
    features_df_marked.loc[:,'end(max)'] = tads_subarray_end
    features_df_marked.loc[:,'chr'] = tads_subarray_chr
    cols = ['chr','start(min)','end(max)']+cols
    features_df_marked = features_df_marked[cols]
    return features_df_marked

def custom_discrete_cmap(numcate,colors,interval):
    newcolors = plt.get_cmap('viridis',numcate).colors
    # assign new colors of interest
    for i in range(0,numcate,interval):
        newcolors[i:i+interval, :] = matplotlib.colors.to_rgba(colors[i])
    # create the customized color map
    cmap = matplotlib.colors.ListedColormap(newcolors)
    return cmap

def plot_heatmap(df,cmap,title,outdir,outname):
    plt.figure()
    ax = sns.heatmap(df,yticklabels=False,cmap=cmap,cbar=False)
    ax.set_title(title)
    ax.vlines([1, 2], colors='black',*ax.get_ylim())
    plt.savefig('{}/{}_{}.png'.format(outdir,outname,title))
    plt.close()

def get_lower_tri_heatmap(df, show_text=False, output="cooc_matrix.png"):
    mask = np.triu(np.ones_like(df))
    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns_plot = sns.heatmap(df, cmap=cmap, vmin=0.5 ,vmax=1, center=0.75, annot=df, mask=mask,
                           square=True, linewidths=.5, cbar_kws={"shrink": .5},annot_kws={"fontsize":16,"color":"black"})
    for text in ax.texts:
        text.set_visible(show_text)
    plt.savefig(output)
    plt.close()

def get_tri_heatmap(df, loc='lower', vmin=0.5, vmax=1, fontsize=16, show_text=False, output="cooc_matrix.png"):
    if loc == 'lower':
        df_tri = df.where(np.tril(np.ones(df.shape)).astype(bool))
    elif loc == 'upper':
        df_tri = df.where(np.triu(np.ones(df.shape)).astype(bool))
    else:
        print("Wrong loc {} params, use lower".format(loc))
        df_tri = df.where(np.tril(np.ones(df.shape)).astype(bool))
    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns_plot = sns.heatmap(df_tri, cmap=cmap, vmin=vmin, vmax=vmax, center=(vmin+vmax)/2, annot=df,
                           square=True, linewidths=.5, cbar_kws={"shrink": .5},
                           annot_kws={"fontsize": fontsize, "color": "black"})
    if loc == 'upper':
        ax.xaxis.tick_top()  # x axis on top
        ax.xaxis.set_label_position('top')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')

    for text in ax.texts:
        text.set_visible(show_text)
    plt.savefig(output)
    plt.close()

def summarize_changes(array_dict,changetype,ctrl_samples,trt_samples):
    summarize_df = pd.DataFrame(0, index=trt_samples, columns=['N','D','S','F'])
    for cs in ctrl_samples:
        cs_array = array_dict[cs]
        for ts in trt_samples:
            ts_array = array_dict[ts]
            if changetype == 'ND':
                if len(ts_array[ts_array<0]) > len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'D'] += 1
                elif len(ts_array[ts_array<0]) < len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'N'] += 1
                else:
                    summarize_df.loc[ts, 'N'] += 0
            elif changetype == 'SF':
                if len(ts_array) > len(cs_array):
                    summarize_df.loc[ts, 'S'] += 1
                elif len(ts_array) < len(cs_array):
                    summarize_df.loc[ts, 'F'] += 1
                else:
                    summarize_df.loc[ts, 'F'] += 0
            elif changetype == 'Mixed':
                if len(ts_array[ts_array<0]) > len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'D'] += 1
                elif len(ts_array[ts_array<0]) < len(cs_array[cs_array<0]):
                    summarize_df.loc[ts, 'N'] += 1
                else:
                    summarize_df.loc[ts, 'N'] += 0
                if len(ts_array) > len(cs_array):
                    summarize_df.loc[ts, 'S'] += 1
                elif len(ts_array) < len(cs_array):
                    summarize_df.loc[ts, 'F'] += 1
                else:
                    summarize_df.loc[ts, 'F'] += 0
    return summarize_df

from functools import reduce
import pandas as pd

def summarize_change_all(dchange_df, tads_sub_array_dict, ctrl_samples, trt_samples):
    summarize_list_dict = {'HS': {'ND': [], 'SF': [], 'Mixed': []},
                           'LS': {'ND': [], 'SF': [], 'Mixed': []}}
    summarize_dict = {}
    sample_changes = {sample: [] for sample in trt_samples}
    sample_changes['tads_id'] = []

    for idx, row in dchange_df.iterrows():
        current_summarize = summarize_changes(tads_sub_array_dict[row['tads_id']], row['ChangeTypes'], ctrl_samples, trt_samples)
        sample_changes['tads_id'].append(row['tads_id'])
        for sample in trt_samples:
            sample_changes[sample].append("{}N;{}D;{}S;{}F".format(
                current_summarize.loc[sample, 'N'],
                current_summarize.loc[sample, 'D'],
                current_summarize.loc[sample, 'S'],
                current_summarize.loc[sample, 'F']
            ))

        summarize_list_dict[row['IndividualChange']][row['ChangeTypes']].append(current_summarize)

        changes = '{}N;{}D;{}S;{}F'.format(
            current_summarize['N'].sum(),
            current_summarize['D'].sum(),
            current_summarize['S'].sum(),
            current_summarize['F'].sum()
        )
        summarize_dict[row['tads_id']] = changes

    def reduce_if_not_empty(lst):
        return reduce(lambda x, y: x.add(y, fill_value=0), lst) if lst else pd.DataFrame()

    summarize_all_dict = {
        'HS': {
            'ND': reduce_if_not_empty(summarize_list_dict['HS']['ND']),
            'SF': reduce_if_not_empty(summarize_list_dict['HS']['SF']),
            'Mixed': reduce_if_not_empty(summarize_list_dict['HS']['Mixed'])
        },
        'LS': {
            'ND': reduce_if_not_empty(summarize_list_dict['LS']['ND']),
            'SF': reduce_if_not_empty(summarize_list_dict['LS']['SF']),
            'Mixed': reduce_if_not_empty(summarize_list_dict['LS']['Mixed'])
        }
    }

    return summarize_dict, summarize_all_dict, pd.DataFrame(sample_changes)


def plotlogo(summarize_df,outname):
    # visualize summarize_df
    color_scheme = {
        'N': [0, .5, 0],
        'D': [1, 0, 0],
        'S': [1, .65, 0],
        'F': [0, 0, 1]
    }
    tick_labels = summarize_df.index.values
    summarize_df_norm = summarize_df/summarize_df.sum(axis=1).max()
    summarize_df_norm.reset_index(drop=True,inplace=True)
    logo = logomaker.Logo(summarize_df_norm, font_name='Arial Rounded MT Bold',color_scheme=color_scheme)
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left'], visible=True, bounds=[0, 2])
    logo.style_spines(spines=['bottom'], visible=True)
    logo.ax.set_ylabel(r"Scaled"+"\n"+ r"Frequency",fontsize=12,weight='bold')
    logo.ax.set_ylim([0 , 1])
    logo.ax.set_yticks([0, 1])
    logo.ax.set_xticks([])
    tick_locations = range(len(tick_labels))
    logo.ax.set_xticks(tick_locations)
    logo.ax.set_xticklabels(tick_labels)
    logo.ax.tick_params(axis='both', which='major', labelsize=16, labelcolor='black')
    for label in logo.ax.get_xticklabels() + logo.ax.get_yticklabels():
        label.set_weight("bold")
    plt.savefig(outname,dpi=500)
    plt.close()
    return summarize_df_norm

def plotHeatmap(features_df_marked, types):
    # individual specific number
    # indi_spe_list = ['high_high_high','high_medium_high','high_low_high']
    # plot clustermap
    value2int = {'high':3.0,'medium':2.0,'low':1.0}
    # high: palevioletred, low:tan, con:forestgreen
    features_list = []
    features_len_list = []
    for group,df in features_df_marked.groupby('GroupChange'):
        print(group)
        if group in ['C','MV']:
            current_df = df[types].replace(value2int)
            current_df.sort_values([types],inplace=True)
            features_list.append(current_df)
            features_list.append(pd.DataFrame(0,index=np.arange(30),columns=list(current_df.columns)))
            features_len_list.append(len(current_df))
        else:
            current_indi = df[(df[types[0]]=='high')&(df[types[2]]=='high')][[types]].replace(value2int)
            current_indi.sort_values(types,inplace=True)
            features_list.append(current_indi)
            features_list.append(pd.DataFrame(0,index=np.arange(30),columns=list(current_indi.columns)))
            features_len_list.append(len(current_indi))
            current_noindi = df[(df[types[0]]=='high')&(df[types[2]]!='high')][types].replace(value2int)
            current_noindi.sort_values(types,inplace=True)
            features_list.append(current_noindi)
            features_len_list.append(len(current_noindi))

    features_final_df = pd.concat(features_list).T
    cmap = custom_discrete_cmap(4,['white','limegreen','goldenrod','red'],1)
    plt.figure(figsize=(15,6))
    ax = sns.heatmap(features_final_df, xticklabels=False, cmap=cmap)
    ax.set_title(types[0])
    ax.hlines([1, 2], colors='black', *ax.get_xlim())
    plt.savefig(f'{outdir}/{types[0]}_heatmap.png')
    plt.close()

def expand_dict_elements(input_dict, rep=5):
    expanded_dict = {}
    for key, value in input_dict.items():
        expanded_dict[key] = [
            [f"{element}{i}" for element in sublist for i in range(rep)]
            for sublist in value
        ]
    return expanded_dict

