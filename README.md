
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dncitPaper

<!-- badges: start -->
<!-- badges: end -->

The goal of dncitPaper is to reproduce the code in the paper Deep
Nonparametric Conditional Independence Tests. To facilitate this
process, this repository is also available as an R package. Note that
access to the UK Biobank (UKB) brain imaging data and its associated
metadata is required to reproduce the study results. Additionally, the
R-package `DNCIT` has to be installed.

## Installation

You can install the development version of dncitPaper from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MSimnach/dncitPaper")
```

## Reproducing the Results

The paper includes three empirical studies:

1.  **Simulation Study**: Located in the folder `inst/simulation_study`
2.  **Study on Brain Structures and Behavioral Traits**: Implemented in
    the files `inst/brain_and_behavior/ukb_preprocess.R` and
    `inst/brain_and_behavior/brain_and_behavior.R`
3.  **Confounder Control Study**: Found in the files
    `inst/brain_and_behavior/ukb_preprocess.R` and
    `inst/brain_and_behavior/confounding_brain_and_behavior.R`.

# Data

To replicate the results, the variable IDs for the UK Biobank data are
available in the folder `inst/extdata/ids`.

# Simulation Study

For the simulation study, datasets need to be generated containing the
IDs of the subjects and the feature representations $X^\omega$ from
Fastsurfer, Freesurfer, and the conditional VAE embedding map.
Additionally, datasets with the IDs of subjects and the confounders (one
dataset per set of confounders) should be prepared. Subjects with
missing values in any of these variables should be excluded. The script
for this process is located in
`inst/simulation_study/sim_XZ_data_gen.R`.

The simulation study can be executed in `R` using the script
`inst/simulation_study/sim_run_settings.R`. Within this script, you can
configure settings such as the DNCIT, the embedding map, the
confounders, their relationship to $Y$, and conditional independence.
Simulations are then run in parallel using the file
`inst/simulation_study/sim_ukb_brainmri.R`, and the results (p-values,
rejection rates, and runtime) are automatically saved after each
configuration.

**Before running the simulation study**, ensure that the path to the UKB
data is correctly specified in the file `R/data_gen.R` at the start of
the function `data_gen()`. The datasets generated at the start should be
saved as in the functions `load_X_orig`, `load_X_obs` and `load_Z`.

# Real-world Applications

In both real-world applications, the UKB data has to be preprocessed
first using the script `inst/brain_and_behavior/ukb_preprocess.R`.
