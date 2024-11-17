#!/usr/bin/env python
#--coding:utf-8 --


#################################################################################################
#################################################################################################
########                                                                                 ########
########    Group-, Individual-specific sample TADs analysis                             ########
########                                                                                 ########
########    Author:  Kun Fang                                                            ########
########                                                                                 ########
########                                                                                 ########
########    Working Environment:  Python3                                                ########
########                                                                                 ########
########    Date:      2024-11-13                                                        ########
########                                                                                 ########
########                                                                                 ########
#################################################################################################
#################################################################################################

import click
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from .GISTA_Util import *

# class for setting up subcommand order.
# Inspired by https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
class SpecialHelpOrder(click.Group):

    def __init__(self, *args, **kwargs):
        self.help_priorities = {}
        super(SpecialHelpOrder, self).__init__(*args, **kwargs)

    def get_help(self, ctx):
        self.list_commands = self.list_commands_for_help
        return super(SpecialHelpOrder, self).get_help(ctx)

    def list_commands_for_help(self, ctx):
        """reorder the list of commands when listing the help"""
        commands = super(SpecialHelpOrder, self).list_commands(ctx)
        return (c[1] for c in sorted(
            (self.help_priorities.get(command, 1), command)
            for command in commands))

    def command(self, *args, **kwargs):
        """Behaves the same as `click.Group.command()` except capture
        a priority for listing command names in help.
        """
        help_priority = kwargs.pop('help_priority', 1)
        help_priorities = self.help_priorities

        def decorator(f):
            cmd = super(SpecialHelpOrder, self).command(*args, **kwargs)(f)
            help_priorities[cmd.name] = help_priority
            return cmd

        return decorator

class Config(object):

    def __init__(self):
        self.verbose = False

pass_config = click.make_pass_decorator(Config, ensure=True)
@click.group(cls=SpecialHelpOrder)
@click.version_option(version='0.1')
@pass_config
def cli(config):
    pass

error_message = "Please contact author kf2799@cumc.columbia.edu for further supporting, thanks for your interesting!"

@cli.command(help_priority=1, help='Multi-samples mode for GISTA')
@click.option('--samplesfile', '-sf', type=click.Path(exists=True),
              help='sample metadata sheet, contain TopDom (like) file list.')
@click.option('--comparison', '-c', type=str, required=True, help="Comparison string in the format 'RT,PT;PT-RT,NT', "
                                                                  "The comparison separate by ';' and treat and control separate by ',',"
                                                                  "The second sample is control")
@click.option('--binsize', '-bs', default=40000, type=int, help = 'resolution/binsize of the TADs')
@click.option('--groupcut', '-gc', default=0.7, type=float, help = 'Group level high variation cutoff')
@click.option('--individualcut', '-ic', default=0.7, type=float, help = 'Individual level high variation cutoff')
@click.option('--outdir', '-od', type=click.Path(), help = 'Output folder')
def multi(samplesfile, comparison, binsize, group_cutoff, individual_cutoff, outdir):
    if outdir is None:
        outdir = './'
    samples_df = pd.read_excel(samplesfile)
    comparison_dict, comparison_types = generate_comparison_dict(samples_df, comparison)
    samples = samples_df['IndividualID'].values
    tads_files = samples_df['FileName'].values
    Path(outdir).mkdir(parents=True, exist_ok=True)

    Tads_array = TadsArray(tads_files,samples,binsize)
    data_tads, tad_ids = Tads_array.loading_files()
    data_dist_tads_df,data_stat_tads_df = Tads_array.tads_basic_stats(data_tads)
    data_stat_tads_df.to_csv('{}/tads_basic_information.csv'.format(outdir),index=False)
    # the format of the key in tads_array_dict: 'chr1_1' means the first tad array in chr1
    # key in sample_df_dict: chr1_1; value: chr/start/end/detail_info
    tads_sub_array_dict, tads_sub_annot_dict, tads_sub_id_list = Tads_array.tads_array_chunk(data_tads,tad_ids)
    # tads_array_size
    tads_sub_array_size = []
    for tads_id in tads_sub_id_list:
        current_size = 0
        for sample in samples:
            current_size += (tads_sub_annot_dict[tads_id][sample][2] - tads_sub_annot_dict[tads_id][sample][1])/1000000
        current_size /= len(samples)
        tads_sub_array_size.append(current_size)

    plt.figure()
    g = sns.displot(data=pd.DataFrame(tads_sub_array_size), x=0, log_scale=10)
    plt.xlabel('TADs array size (Mb)')
    plt.savefig('{}/TADs_sub_array_size_dist.png'.format(outdir))
    plt.tight_layout()
    plt.close()

    Tads_matrix = TadsMatrix(tads_sub_array_dict,samples,tads_sub_id_list)
    # trt type at 0 position and ctrl type at the 1 position
    # total_info with key total_detail, total_sf, total_nd, nd:Neo-Del, sf: Split-Fuse
    total_info, features_df = Tads_matrix.fit(comparison_dict)
    features_df.to_csv("{}/Tads_sub_array_feature_matrix.csv".format(outdir),index=False)
    features_dict = {}
    for comp in comparison_types:
        cur_types = comparison_types[comp]
        nd_type = [comp_type + '_nd' for comp_type in cur_types]
        sf_type = [comp_type + '_sf' for comp_type in cur_types]
        features_nd = features_df.loc[:, nd_type]
        features_sf = features_df.loc[:, sf_type]
        # calculate pval for SV
        print("Permutation Test")
        cur_pvals = cal_pval(pd.concat([features_nd, features_sf], axis=1), comp)
        # give mnemonic
        features_nd_marked, nd_cuts = class_process(features_nd, nd_type[0], nd_type[2], nd_type[1], group_cutoff, individual_cutoff)
        features_sf_marked, sf_cuts = class_process(features_sf, sf_type[0], sf_type[2], sf_type[1], group_cutoff, individual_cutoff)
        # nd sf overall distribution
        nd_data = features_nd.stack()
        sf_data = features_sf.stack()
        # visual 
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        f, ax = plt.subplots(1, 2, sharey=True, sharex=True)
        ax1 = ax[0]
        sns.kdeplot(data=np.log2(nd_data + 0.1), ax=ax1, color='black', fill=True)
        ax1.axvline(np.log2(nd_cuts[0] + 0.1), linewidth=1, color='r', ls='--')
        ax1.axvline(np.log2(nd_cuts[1] + 0.1), linewidth=1, color='g', ls='--')
        ax1.set_ylabel(f'{comp}\nDensity')
        ax1.set_xlabel('Log2(Neo-Del Score+0.1)')
        ax2 = ax[1]
        sns.kdeplot(data=np.log2(sf_data + 0.1), ax=ax2, color='purple', fill=True)
        ax2.axvline(np.log2(sf_cuts[0] + 0.1), linewidth=1, color='r', ls='--')
        ax2.axvline(np.log2(sf_cuts[1] + 0.1), linewidth=1, color='g', ls='--')
        ax2.set_xlabel('Log2(Split-Fuse Score+0.1)')
        plt.tight_layout()
        plt.savefig(f'{outdir}/ND_SF_toCate.png', dpi=300)
        plt.close()
        
        # get annotation
        features_df_marked = combine_nd_sf_annot(features_df['tads_id'], features_nd_marked,
                                                        features_sf_marked,
                                                        cur_types[0], cur_types[2], cur_types[1])
        # assign coordinate back
        features_df_marked = coor_annot(features_df_marked, tads_sub_id_list, tads_sub_annot_dict)
        # differential changes, specific change type
        features_dchange = features_df_marked.loc[features_df_marked['GroupChange'] == 'SV', :]
        summarize_dict, summarize_list_dict, sample_changes = summarize_change_all(
            features_dchange,
            tads_sub_array_dict, comparison_dict[comparison_types[comp][2]][1],
            comparison_dict[comparison_types[comp][1]][1])
        # the maximum changecode's number should be less than len(ctrl_samples)*len(trt_samples)
        features_df_marked['ChangeCodes'] = features_df_marked['tads_id'].map(summarize_dict)
        features_df_marked = pd.merge(features_df_marked, sample_changes, on='tads_id',how='outer')
        features_df_marked = features_df_marked.fillna('NotAvail')
        # plot sample change logo
        ivh_nd_norm = plotlogo(summarize_list_dict['HS']['ND'], '{}/{}_IVH_ND_logo.png'.format(outdir, comp))
        ivh_sf_norm = plotlogo(summarize_list_dict['HS']['SF'], '{}/{}_IVH_SF_logo.png'.format(outdir, comp))
        ivh_mix_norm = plotlogo(summarize_list_dict['HS']['Mixed'], '{}/{}_IVH_Mixed_logo.png'.format(outdir, comp))

        imc_nd_norm = plotlogo(summarize_list_dict['LS']['ND'], '{}/{}_IMC_ND_logo.png'.format(outdir, comp))
        imc_sf_norm = plotlogo(summarize_list_dict['LS']['SF'], '{}/{}_IMC_SF_logo.png'.format(outdir, comp))
        imc_mix_norm = plotlogo(summarize_list_dict['LS']['Mixed'], '{}/{}_IMC_Mixed_logo.png'.format(outdir, comp))

        nd_norm = plotlogo(summarize_list_dict['HS']['ND'] + summarize_list_dict['LS']['ND'],
                           '{}/{}_ND_logo.png'.format(outdir, comp))
        sf_norm = plotlogo(summarize_list_dict['HS']['SF'] + summarize_list_dict['LS']['SF'],
                           '{}/{}_SF_logo.png'.format(outdir, comp))
        mix_norm = plotlogo(summarize_list_dict['HS']['Mixed'] + summarize_list_dict['LS']['Mixed'],
                            '{}/{}_Mixed_logo.png'.format(outdir, comp))

        # finalize the pval
        mask = np.full(len(features_df_marked), False)
        mask[features_df_marked[features_df_marked['GroupChange'] == "SV"].index.values] = True
        features_df_marked['pval_SV'] = np.where(mask, np.minimum(cur_pvals[0], cur_pvals[1]),cur_pvals[1])
        features_df_dchange = features_df_marked.loc[features_df_marked['GroupChange'] == 'SV', :]
        # save as excel
        with pd.ExcelWriter("{}/{}_tads_sub_array_marks.xlsx".format(outdir,comp),
                            engine="xlsxwriter") as writer:
            features_df_marked.to_excel(writer, index=False, sheet_name=f"{comp}_total_info")
            features_df_dchange.to_excel(writer, index=False, sheet_name=f"{comp}_SV")

        features_TTvsNT_df_marked['GroupChange'] = pd.Categorical(features_TTvsNT_df_marked['GroupChange'],
                                                                  categories=['C', 'MV', 'SV'], ordered=True)
        plotHeatmap(features_df_marked, [cur_types[0], cur_types[2], cur_types[1]])
        features_dict[comp] = features_df_marked
    
    # summary
    summary_dict = {'GroupChange': [], 'Type': [], 'Count': []}
    for comp in comparison_dict:
        cur_feature_df_marked = features_dict[comp]
        bio_mark = Counter(cur_feature_df_marked['GroupChange'])
        for key in bio_mark:
            summary_dict['GroupChange'].append(key)
            summary_dict['Type'].append(comp)
            summary_dict['Count'].append(bio_mark[key])

    summary_df = pd.DataFrame(summary_dict)
    plt.figure()
    sns.barplot(data=summary_df,x='GroupChange',y='Count',hue='Type',order=['C','MV','SV'])
    plt.savefig('{}/TADs_array_GroupChange_summary.png'.format(outdir))
    plt.close()

    # calculate correlations for all sample based on tads_array
    total_size = 0
    tadarray_size_dict = {}
    for tadarray in tads_sub_array_dict:
        current_tadarray = tads_sub_array_dict[tadarray]
        max_size = max((sum(np.abs(v))) for k, v in current_tadarray.items())
        total_size += max_size
        tadarray_size_dict[tadarray] = max_size
    sample_tad_corr_dict = {}
    for tadarray in tads_sub_array_dict:
        current_tadarray = tads_sub_array_dict[tadarray]
        # pearson correlation should have same length and at least 2 length
        max_length = max((len(v)) for k, v in current_tadarray.items())
        max_size = tadarray_size_dict[tadarray]
        tadarray_weight = max_size / total_size
        if max_length == 1:
            max_length = 2
        new_tadarray = {}
        for sample in current_tadarray:
            sample_array = current_tadarray[sample].tolist() + (max_length - len(current_tadarray[sample])) * [0]
            # avoid NAN generated based no variance within sample, e.g., [1,1,1]
            if len(np.unique(sample_array)) == 1:
                sample_array[0] -= 0.1
                sample_array[1] += 0.1
            new_tadarray[sample] = sample_array
        samples_tad_corr_mat = pd.DataFrame(0, columns=samples, index=samples)
        for subset in itertools.combinations(samples, 2):
            sample1 = subset[0]
            sample2 = subset[1]
            pearsonr_result = pearsonr(new_tadarray[sample1],
                                       new_tadarray[sample2])
            samples_tad_corr_mat.loc[sample1, sample2] = pearsonr_result[0]
            samples_tad_corr_mat.loc[sample2, sample1] = pearsonr_result[0]
        # > 0.7 to 1; >0.3 <0.7 to 0.5; <0.3 to 0
        samples_tad_corr_mat[samples_tad_corr_mat >= 0.3] = 1
        samples_tad_corr_mat[(samples_tad_corr_mat > -0.3) & (samples_tad_corr_mat < 0.3)] = 0.5
        samples_tad_corr_mat[samples_tad_corr_mat <= -0.3] = 0
        sample_tad_corr_dict[tadarray] = tadarray_weight * samples_tad_corr_mat

    samples_tad_corr_mat_final = sum(sample_tad_corr_dict.values())
    samples_tad_corr_mat_final.values[np.diag_indices_from(samples_tad_corr_mat_final)] = 1
    get_tri_heatmap(samples_tad_corr_mat_final, loc='lower', show_text=False, vmin=0.6, vmax=0.9,
                    output='{}/samples_pearsonr_heatmap.png'.format(outdir))


@cli.command(help_priority=1, help='Multi-samples mode for GISTA')
@click.option('--samplesfile', '-sf', type=click.Path(exists=True),
              help='sample metadata sheet, contain TopDom (like) file list.')
@click.option('--comparison', '-c', type=str, required=True, help="Comparison string in the format 'RT,PT;PT-RT,NT', "
                                                                  "The comparison separate by ';' and treat and control separate by ',',"
                                                                  "The second sample is control")
@click.option('--binsize', '-bs', default=40000, type=int, help = 'resolution/binsize of the TADs')
@click.option('--groupcut', '-gc', default=0.3, type=float, help = 'Group level high variation cutoff')
@click.option('--individualcut', '-ic', default=0.7, type=float, help = 'Individual level high variation cutoff')
@click.option('--pseudorep', '-pr', default=5, type=int, help = 'The number of Pseudo-replication')
@click.option('--outdir', '-od', type=click.Path(), help = 'Output folder')
def two(samplesfile, comparison, binsize, group_cutoff, individual_cutoff, pseudorep, outdir):
    if outdir is None:
        outdir = './'
    samples_df = pd.read_excel(samplesfile)
    samples = samples_df['IndividualID'].values
    tads_files = samples_df['FileName'].values
    Path(outdir).mkdir(parents=True, exist_ok=True)

    Tads_array = TadsArray(tads_files, samples, binsize)
    data_tads, tad_ids = Tads_array.loading_files()
    data_dist_tads_df, data_stat_tads_df = Tads_array.tads_basic_stats(data_tads)
    data_stat_tads_df.to_csv('{}/tads_basic_information.csv'.format(outdir), index=False)
    # the format of the key in tads_array_dict: 'chr1_1' means the first tad array in chr1
    # key in sample_df_dict: chr1_1; value: chr/start/end/detail_info
    tads_sub_array_dict, tads_sub_annot_dict, tads_sub_id_list = Tads_array.tads_array_chunk(data_tads, tad_ids)
    # tads_array_size
    tads_sub_array_size = []
    for tads_id in tads_sub_id_list:
        current_size = 0
        for sample in samples:
            current_size += (tads_sub_annot_dict[tads_id][sample][2] - tads_sub_annot_dict[tads_id][sample][
                1]) / 1000000
        current_size /= len(samples)
        tads_sub_array_size.append(current_size)

    plt.figure()
    g = sns.displot(data=pd.DataFrame(tads_sub_array_size), x=0, log_scale=10)
    plt.xlabel('TADs array size (Mb)')
    plt.savefig('{}/TADs_sub_array_size_dist.png'.format(outdir))
    plt.tight_layout()
    plt.close()

    Tads_matrix = TadsMatrix(tads_sub_array_dict, samples, tads_sub_id_list,onerep=True, pseudorep=pseudorep)
    comparison_tmp_dict, comparison_types = generate_comparison_dict(samples_df, comparison)
    comparison_dict = expand_dict_elements(comparison_tmp_dict, pseudorep)
    # trt type at 0 position and ctrl type at the 1 position
    # total_info with key total_detail, total_sf, total_nd, nd:Neo-Del, sf: Split-Fuse
    total_info, features_df = Tads_matrix.fit(comparison_dict)
    features_df.to_csv("{}/Tads_sub_array_feature_matrix.csv".format(outdir), index=False)
    features_dict = {}
    for comp in comparison_types:
        cur_types = comparison_types[comp]
        nd_type = [comp_type + '_nd' for comp_type in cur_types]
        sf_type = [comp_type + '_sf' for comp_type in cur_types]
        features_nd = features_df.loc[:, nd_type]
        features_sf = features_df.loc[:, sf_type]
        # calculate pval for SV
        print("Permutation Test")
        cur_pvals = cal_pval(pd.concat([features_nd, features_sf], axis=1), comp)
        # give mnemonic
        features_nd_marked, nd_cuts = class_process(features_nd, nd_type[0], nd_type[2], nd_type[1], group_cutoff,
                                                    individual_cutoff)
        features_sf_marked, sf_cuts = class_process(features_sf, sf_type[0], sf_type[2], sf_type[1], group_cutoff,
                                                    individual_cutoff)
        # nd sf overall distribution
        nd_data = features_nd.stack()
        sf_data = features_sf.stack()
        # visual 
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        f, ax = plt.subplots(1, 2, sharey=True, sharex=True)
        ax1 = ax[0]
        sns.kdeplot(data=np.log2(nd_data + 0.1), ax=ax1, color='black', fill=True)
        ax1.axvline(np.log2(nd_cuts[0] + 0.1), linewidth=1, color='r', ls='--')
        ax1.axvline(np.log2(nd_cuts[1] + 0.1), linewidth=1, color='g', ls='--')
        ax1.set_ylabel(f'{comp}\nDensity')
        ax1.set_xlabel('Log2(Neo-Del Score+0.1)')
        ax2 = ax[1]
        sns.kdeplot(data=np.log2(sf_data + 0.1), ax=ax2, color='purple', fill=True)
        ax2.axvline(np.log2(sf_cuts[0] + 0.1), linewidth=1, color='r', ls='--')
        ax2.axvline(np.log2(sf_cuts[1] + 0.1), linewidth=1, color='g', ls='--')
        ax2.set_xlabel('Log2(Split-Fuse Score+0.1)')
        plt.tight_layout()
        plt.savefig(f'{outdir}/ND_SF_toCate.png', dpi=300)
        plt.close()

        # get annotation
        features_df_marked = combine_nd_sf_annot(features_df['tads_id'], features_nd_marked,
                                                 features_sf_marked,
                                                 cur_types[0], cur_types[2], cur_types[1])
        # assign coordinate back
        features_df_marked = coor_annot(features_df_marked, tads_sub_id_list, tads_sub_annot_dict)
        # differential changes, specific change type
        features_dchange = features_df_marked.loc[features_df_marked['GroupChange'] == 'SV', :]
        summarize_dict, summarize_list_dict, sample_changes = summarize_change_all(
            features_dchange,
            tads_sub_array_dict, comparison_dict[comparison_types[comp][2]][1],
            comparison_dict[comparison_types[comp][1]][1])
        # the maximum changecode's number should be less than len(ctrl_samples)*len(trt_samples)
        features_df_marked['ChangeCodes'] = features_df_marked['tads_id'].map(summarize_dict)
        features_df_marked = pd.merge(features_df_marked, sample_changes, on='tads_id', how='outer')
        features_df_marked = features_df_marked.fillna('NotAvail')
        # plot sample change logo
        ivh_nd_norm = plotlogo(summarize_list_dict['HS']['ND'], '{}/{}_IVH_ND_logo.png'.format(outdir, comp))
        ivh_sf_norm = plotlogo(summarize_list_dict['HS']['SF'], '{}/{}_IVH_SF_logo.png'.format(outdir, comp))
        ivh_mix_norm = plotlogo(summarize_list_dict['HS']['Mixed'], '{}/{}_IVH_Mixed_logo.png'.format(outdir, comp))

        imc_nd_norm = plotlogo(summarize_list_dict['LS']['ND'], '{}/{}_IMC_ND_logo.png'.format(outdir, comp))
        imc_sf_norm = plotlogo(summarize_list_dict['LS']['SF'], '{}/{}_IMC_SF_logo.png'.format(outdir, comp))
        imc_mix_norm = plotlogo(summarize_list_dict['LS']['Mixed'], '{}/{}_IMC_Mixed_logo.png'.format(outdir, comp))

        nd_norm = plotlogo(summarize_list_dict['HS']['ND'] + summarize_list_dict['LS']['ND'],
                           '{}/{}_ND_logo.png'.format(outdir, comp))
        sf_norm = plotlogo(summarize_list_dict['HS']['SF'] + summarize_list_dict['LS']['SF'],
                           '{}/{}_SF_logo.png'.format(outdir, comp))
        mix_norm = plotlogo(summarize_list_dict['HS']['Mixed'] + summarize_list_dict['LS']['Mixed'],
                            '{}/{}_Mixed_logo.png'.format(outdir, comp))

        # finalize the pval
        mask = np.full(len(features_df_marked), False)
        mask[features_df_marked[features_df_marked['GroupChange'] == "SV"].index.values] = True
        features_df_marked['pval_SV'] = np.where(mask, np.minimum(cur_pvals[0], cur_pvals[1]), cur_pvals[1])
        features_df_dchange = features_df_marked.loc[features_df_marked['GroupChange'] == 'SV', :]
        # save as excel
        with pd.ExcelWriter("{}/{}_tads_sub_array_marks.xlsx".format(outdir, comp),
                            engine="xlsxwriter") as writer:
            features_df_marked.to_excel(writer, index=False, sheet_name=f"{comp}_total_info")
            features_df_dchange.to_excel(writer, index=False, sheet_name=f"{comp}_SV")

        features_TTvsNT_df_marked['GroupChange'] = pd.Categorical(features_TTvsNT_df_marked['GroupChange'],
                                                                  categories=['C', 'MV', 'SV'], ordered=True)
        plotHeatmap(features_df_marked, [cur_types[0], cur_types[2], cur_types[1]])
        features_dict[comp] = features_df_marked

    # summary
    summary_dict = {'GroupChange': [], 'Type': [], 'Count': []}
    for comp in comparison_dict:
        cur_feature_df_marked = features_dict[comp]
        bio_mark = Counter(cur_feature_df_marked['GroupChange'])
        for key in bio_mark:
            summary_dict['GroupChange'].append(key)
            summary_dict['Type'].append(comp)
            summary_dict['Count'].append(bio_mark[key])

    summary_df = pd.DataFrame(summary_dict)
    plt.figure()
    sns.barplot(data=summary_df, x='GroupChange', y='Count', hue='Type', order=['C', 'MV', 'SV'])
    plt.savefig('{}/TADs_array_GroupChange_summary.png'.format(outdir))
    plt.close()

    # calculate correlations for all sample based on tads_array
    total_size = 0
    tadarray_size_dict = {}
    for tadarray in tads_sub_array_dict:
        current_tadarray = tads_sub_array_dict[tadarray]
        max_size = max((sum(np.abs(v))) for k, v in current_tadarray.items())
        total_size += max_size
        tadarray_size_dict[tadarray] = max_size
    sample_tad_corr_dict = {}
    for tadarray in tads_sub_array_dict:
        current_tadarray = tads_sub_array_dict[tadarray]
        # pearson correlation should have same length and at least 2 length
        max_length = max((len(v)) for k, v in current_tadarray.items())
        max_size = tadarray_size_dict[tadarray]
        tadarray_weight = max_size / total_size
        if max_length == 1:
            max_length = 2
        new_tadarray = {}
        for sample in current_tadarray:
            sample_array = current_tadarray[sample].tolist() + (max_length - len(current_tadarray[sample])) * [0]
            # avoid NAN generated based no variance within sample, e.g., [1,1,1]
            if len(np.unique(sample_array)) == 1:
                sample_array[0] -= 0.1
                sample_array[1] += 0.1
            new_tadarray[sample] = sample_array
        samples_tad_corr_mat = pd.DataFrame(0, columns=samples, index=samples)
        for subset in itertools.combinations(samples, 2):
            sample1 = subset[0]
            sample2 = subset[1]
            pearsonr_result = pearsonr(new_tadarray[sample1],
                                       new_tadarray[sample2])
            samples_tad_corr_mat.loc[sample1, sample2] = pearsonr_result[0]
            samples_tad_corr_mat.loc[sample2, sample1] = pearsonr_result[0]
        # > 0.7 to 1; >0.3 <0.7 to 0.5; <0.3 to 0
        samples_tad_corr_mat[samples_tad_corr_mat >= 0.3] = 1
        samples_tad_corr_mat[(samples_tad_corr_mat > -0.3) & (samples_tad_corr_mat < 0.3)] = 0.5
        samples_tad_corr_mat[samples_tad_corr_mat <= -0.3] = 0
        sample_tad_corr_dict[tadarray] = tadarray_weight * samples_tad_corr_mat

    samples_tad_corr_mat_final = sum(sample_tad_corr_dict.values())
    samples_tad_corr_mat_final.values[np.diag_indices_from(samples_tad_corr_mat_final)] = 1
    get_tri_heatmap(samples_tad_corr_mat_final, loc='lower', show_text=False, vmin=0.6, vmax=0.9,
                    output='{}/samples_pearsonr_heatmap.png'.format(outdir))
