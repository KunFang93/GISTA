# GISTA
A light-weight tool for finding Group-, Individual-sample Specific TADs (GISTA)
GISTA is composed of three steps: 1). Partitioning the whole genome into a series of TADs arrays across all samples; 2). Constructing a feature vector for each TADs array; 3). Annotating the type of TAD change to each TADs arrays and statistical inference.

# WorkFlow
<img src="https://github.com/KunFang93/GISTA/blob/main/Workflow.png" width="900">

# Installation
```
Step1: install mamba from https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
Step2: mamba create -n gista python pandas numpy seaborn scipy tqdm logomaker click openpyxl xlsxwriter natsort
Step3: mamba activate gista
Step4: cd <path to the GISTA>
Step5: pip install --editable .
```

# Quick Start
#### Currently GISTA only accept TopDom output formats (please check the exmaple data). If you need help, please contact author
### Multi-Samples Mode
```
# 'RT,PT;PT-RT,NT' mean process RT.vs.PT and PT+RT.vs.NT
# The comparison separate by ';' and treat and control separate by ',',The second sample is control
# see SampleList_Multi.xlsx exmaple in data
GISTA multi -sf SampleList_Multi.xlsx -c 'RT,PT;PT-RT,NT' --binsize 40000
```

### Two-Samples Mode
```
# 'RT,PT;PT-RT,NT' mean process RT.vs.PT and PT+RT.vs.NT
# The comparison separate by ';' and treat and control separate by ',',The second sample is control
# see SampleList_Multi.xlsx exmaple in data
GISTA two -sf SampleList_Two.xlsx -c 'RT,PT;PT-RT,NT' --binsize 40000
```

# Manual
Please use 'GISTA multi --help' or 'GISTA two --help' to check detailed instruction
```
Welcome to use GISTA :)
Usage: GISTA [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  multi  Multi-samples mode for GISTA
  two    Two-samples mode for GISTA
```
**Multi**
```
Usage: GISTA multi [OPTIONS]

  Multi-samples mode for GISTA

Options:
  -sf, --samplesfile PATH     sample metadata sheet, contain TopDom (like)
                              files, .csv or .xlsx  [required]
  -c, --comparison TEXT       Comparison string in the format 'RT,PT;PT-
                              RT,NT', The comparison separate by ';' and treat
                              and control separate by ',',The second sample is
                              control  [required]
  -bs, --binsize INTEGER      resolution/binsize of the TADs
  -gc, --groupcut FLOAT       Group level high variation cutoff
  -ic, --individualcut FLOAT  Individual level high variation cutoff
  -od, --outdir PATH          Output folder
  --help                      Show this message and exit.
```
**Two**
```
Usage: GISTA two [OPTIONS]

  Two-samples mode for GISTA

Options:
  -sf, --samplesfile PATH     sample metadata sheet, contain TopDom (like)
                              file, .csv or .xlsx  [required]
  -c, --comparison TEXT       Comparison string in the format 'RT,PT;PT-
                              RT,NT', The comparison separate by ';' and treat
                              and control separate by ',',The second sample is
                              control  [required]
  -bs, --binsize INTEGER      resolution/binsize of the TADs
  -gc, --groupcut FLOAT       Group level high variation cutoff
  -ic, --individualcut FLOAT  Individual level high variation cutoff
  -pr, --pseudorep INTEGER    The number of Pseudo-replication
  -od, --outdir PATH          Output folder
  --help                      Show this message and exit.
```

# Citation
If you find this tool useful, please consider cite
```
Choppavarapu, L., Fang, K., Liu, T., & Jin, V. X. (2024). Hi-C profiling in tissues reveals 3D chromatin-regulated breast tumor heterogeneity and tumor-specific looping-mediated biological pathways. bioRxiv, 2024-03.
```
