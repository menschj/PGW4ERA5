# Repository pgw-python

Software to modify ERA5 files to impose a large-scale climate change signal in ERA5 files 
as described e.g. here https://iopscience.iop.org/article/10.1088/1748-9326/ab4438 or also
here https://doi.org/10.1175/JCLI-D-18-0431.1

# General Documentation
To modify the ERA5 files, we need a climate change signal obtained from the difference between two GCM climatologies, HIST and SCEN. The climate change signal is thus SCEN-HIST and referred to as climate delta.

The structure of the repository is built as follows:

The top level directory contains the central scripts to preprocess the GCM climate change signal [step_02_preproc_deltas.py](/step_02_preproc_deltas.py) and to modify the ERA5 files [step_03_apply_to_era.py](/step_03_apply_to_era.py) with the climate change signal.
The subdirectory [step_01_extract_deltas](/step_01_extract_deltas/) contains less generic code that can serve as a template to obtain the GCM climatologies HIST, SCEN, as well as the climate delta SCEN-HIST from raw CMIP6 output. The starting point here is the script [extract_climate_delta.sh](/step_01_extract_deltas/extract_climate_delta.sh). This script has to be adjusted depending on the specific use case.

Note that essential usage-oriented information can be found by running `python step_02_preproc_deltas.py --help` and `python step_03_apply_to_era.py --help`.

Note that users that feed the processed ERA5 files into Int2lm (to run COSMO or ICON) should also modify the variable T_CL in their external parameter file, see postproc_cosmo for more information.

# Software Requirements

The software is written in python 3 and requires multiple python modules. The ennvironment-file **environment.yml** can be used to install a conda environment to run the software. More information about what conda is and how it works: https://docs.conda.io/projects/conda/en/latest/user-guide/index.html#

To install the enviroment, just execute `conda env create -f environment.yml` once conda is installed. 

# Workflows Based on Input Data

**Requeriments**

Annual climate deltas (SCEN-HIST) and a historical climatology (HIST) from a global climate model in either daily or monthly steps.
Climate deltas refer to the difference between the fields predicted by the climate model between two different time periods (usually future and present). If climate model data in the CMOR format (e.g. CMIP simulations) will be used to force the PGW simulations there is a practical [documentation](/Documentations/README_CMOR.md) on which variables are needed.
Template scripts to extract CMIP6 data are given in [step_01_extract_deltas](/step_01_extract_deltas/), e.g. [here](/step_01_extract_deltas/extract_climate_delta.sh).

After computing the raw climate deltas on the GCM grid, run the following scripts:
1) Only if using daily climate deltas instead of monthly (note that is not really recommended), smooth deltas in time: `python step_02_preproc_deltas.py smoothing [...]`
2) Regrid deltas to ERA5 grid: `python step_02_preproc_deltas.py regridding [...]`
3) Modify ERA5 files: `python step_03_apply_to_era.py [...]`
4) There may be some additional steps required for a specific limited-area model. For instance in COSMO, the deep soil temperature climatology has to be adjusted in the external parameters file. [postproc_cosmo](/postproc_cosmo/). 
5) After these steps, the limited-area-model-specific routine to convert ERA5 files to model initial and boundary conditions can be run using the modified ERA5 files as input.

# References
To acknowledge this software cite the following article:

Brogli, R., Heim, C., Mensch, J., Sørland, S. L., & Schär, C. (2023). The pseudo-global-warming (PGW) approach: Methodology, software package PGW4ERA5 v1.1, validation and sensitivity analyses. Geoscientific Model Development, preprint. https://doi.org/10.5194/gmd-2022-167 
