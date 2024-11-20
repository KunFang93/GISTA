# GISTA
A light-weight tool for finding Group-, Individual-sample Specific TADs (GISTA)

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

# Citation
If you find this tool useful, please consider cite
```
Choppavarapu, L., Fang, K., Liu, T., & Jin, V. X. (2024). Hi-C profiling in tissues reveals 3D chromatin-regulated breast tumor heterogeneity and tumor-specific looping-mediated biological pathways. bioRxiv, 2024-03.
```
