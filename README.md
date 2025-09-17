Natalie Rhodes et al. Imaging Neuroscience 2025 (https://doi.org/10.1162/imag_a_00527)
--------------------------------------------------------------------------------------------------------------
All data are published on Zenodo (Nottingham - https://zenodo.org/records/15190020) and will be published on Ontario Brain Institute (SickKids)

Requirements:
- Fieldtrip version 20220906 https://www.fieldtriptoolbox.org/download.php
- Matlab version 2022b
- FSL for use of FLIRT to coregister brains to MNI space

Data file structure (pseudo-BIDs format):
- ---> Data
- -----> BIDS
- -------> derivatives
- -------> sub-001
- -------> sub-002
- -------> sub-00n
- ---------> ses-001
- -----------> anat
- -----------> meg

To run analyses:
- 'run_opm_gamma_analyses.m'

which calls in the following functions: 
- 'data_quality_func.m' : checks data quality producing a table of information
- 'opm_gamma_4mm.m' : runs beamformer on 4mm brain grid on gamma-band filtered data, outputting 4mm resolution pseudo-T statistical maps
- 'opm_gamma_visualmask.m': runs beamformer on 1mm grid on gamma-band filtered data across a dilated mask of the visual cortex, ouputing peak virtual electrode
- 'opm_gamma_4mm_alpha': runs the 4mm brain grid on alpha-band filtered data outputting 4mm resolution pseudo-T statistical maps

shell scripts:
- tstat2mni_gamma.sh
- tstat2mni_alpha.sh
transform pseudo-T maps from individual anatomy to MNI space for averaging

- /plotting/average_pseudoT.m: averages the pseudo-T maps in MNI space to produce group average images
- /plotting/plot_peak_locs_MNI.m plots the peak location mean and standard deviation ellipsoids on MNI mesh

