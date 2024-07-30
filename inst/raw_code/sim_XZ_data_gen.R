#### Extraction of the brain MRIs together with confounders age and sex as well as together with
# age, sex and 10 genetic PCs to load them later more easily during the simulation study
library(data.table)
library(dplyr)
library(stringr)

# paths to data repository
path_to_ukb_data <- 'M:/CITs/Application/UKB_data/ukb49727.csv'
path_to_fastsurfer_ids <- 'M:/CITs/Application/UKB_data/ids/Ids_IDPs.csv'
path_to_save_fastsurfer_X <- "M:\\CITs\\Application\\UKB_data\\ukb_fastsurfer.csv"
path_to_save_freesurfer_X <- "M:\\CITs\\Application\\UKB_data\\ukb_freesurfer.csv"
path_to_save_age_sex_Z <- "M:\\CITs\\Application\\UKB_data\\ukb_Z_age_sex.csv"
path_to_save_age_sex_10_genes_Z <- "M:\\CITs\\Application\\UKB_data\\ukb_Z_genes10.csv"
path_to_confounder_ids <- 'M:/CITs/Application/UKB_data/ids/ids_confounder_avinun.csv'

### Subset available UKB to obtain available col names
ukb_whole_columns <- data.table::fread(file=path_to_ukb_data, header=TRUE, nrows=10)

# check ids from UKB showcase
cols_to_select <- colnames(ukb_whole_columns)[grepl('41202', colnames(ukb_whole_columns), fixed = TRUE)]
ukb_whole_columns[['53-2.0']][9]
id_eid = 'eid'
id_sex <- '31-0.0'
id_age_assessment_center <- '21003-2.0'
id_assessment_center <- '54-0.0'
id_home_location_urban <- '20118-0.0'
ids_genetic_pcs <- colnames(ukb_whole_columns)[grepl('22009', colnames(ukb_whole_columns), fixed = TRUE)]

# Fastsurfer feature representation ids
ids_IDPs <- data.table::fread(file =path_to_fastsurfer_ids, select='Field ID')
ids_IDPs_full_str <- ids_IDPs %>%
  dplyr::mutate_all(list(~ stringr::str_c(., '-2.0'))) %>%
  dplyr::pull(`Field ID`)

ids_ukb_brain_mri <- c(id_eid, id_sex, id_age_assessment_center,id_assessment_center,
                       id_home_location_urban, ids_genetic_pcs,
                       ids_IDPs_full_str)

### load data depending on selected ids
ukb_data <- data.table::fread(file=path_to_ukb_data,select = ids_ukb_brain_mri,  header=TRUE)

ids_confounders <- c(id_eid, id_sex, id_age_assessment_center,id_assessment_center,
                     id_home_location_urban, ids_genetic_pcs)

### Age sex as confounders
ukb_pipeline <- stats::na.omit(ukb_data)
ukb_Z <- ukb_pipeline[,1:3]
data.table::fwrite(ukb_Z, file = path_to_save_age_sex_Z)

### Embedding map: Fastsurfer
n_confounders <- length(ids_confounders)
ukb_fastsurfer <- ukb_pipeline[,c(1, (n_confounders+1):(n_confounders+139)), with=FALSE]
data.table::fwrite(ukb_fastsurfer, file = path_to_save_fastsurfer_X)

### Embedding map: Freesurfer
cols_aseg <- 140:238+n_confounders
cols_ba_exvivo <- 239:322+n_confounders
cols_a2009s <- 323:766+n_confounders
cols_DKT <- 767:952+n_confounders
cols_desikan_gw <- 953:1022+n_confounders
cols_desikan_pial <- 1023:1088+n_confounders
cols_desikan_white <- 1089:1290+n_confounders
cols_subseg <- 1291:1411+n_confounders
cols_freesurfer_selected <- c(cols_aseg, cols_desikan_pial)
ukb_freesurfer <- ukb_pipeline[,c(1, cols_freesurfer_selected), with=FALSE]
data.table::fwrite(ukb_freesurfer, file=path_to_save_freesurfer_X)

### Embedding from conditional VAE
#ukb_transfer <- fread(file='M:/CITs/Application/UKB_data/ukb_condVAE.csv',  header=FALSE)

#### ten dimensional confounders of genetic PCs
ukb_pipeline_ten <- stats::na.omit(ukb_data)
ukb_Z_genes10 <- ukb_pipeline_ten[,c(1:3, 6:15)]
data.table::fwrite(ukb_Z_genes10, file = path_to_save_age_sex_10_genes_Z)


#### Multiple confounder sets
path_to_save_confounders <- "M:\\CITs\\Application\\UKB_data\\"
# get ids of confounders (load additional fastsurfer features to reduce sample size by removing obs. without fastsurfer)
ids_confounders <- data.table::fread(file=path_to_confounder_ids,header=TRUE)$ids_confounder
ids_ukb_brain_mri_confounds <- c(id_eid, ids_confounders,
                                 ids_IDPs_full_str)
ukb_data <- data.table::fread(file=path_to_ukb_data,select = ids_ukb_brain_mri_confounds,  header=TRUE)
# remove obs. with missing values
ukb_wo_missing <- stats::na.omit(ukb_data)
#select only confounders
ukb_confounders <- ukb_wo_missing[,1:16]

## Preprocess confounds (site & date)
colnames(ukb_confounders) <- c('id', 'sex', 'age', 'site', 'date',
                               'head_size', 'head_location1', 'head_location2','head_location3','head_location4',
                               'qc_discrepancy', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5')
## dummy confounder site
dummy_df <- fastDummies::dummy_cols(ukb_confounders, select_columns = 'site')
#remove original site column and one site (used as baseline, else non-identifiable)
remove_sites <- !names(dummy_df) %in% c('site', 'site_10003')
confounds_dummy <- dummy_df[, .SD, .SDcols = remove_sites]
## date difference to first date
confounds_dummy$date_diff <- confounds_dummy[['date']] - min(confounds_dummy[['date']])
# remove original date column
remove_date <- !names(confounds_dummy) %in% c('date')
confounds_dummy <- confounds_dummy[, .SD, .SDcols = remove_date]
# confounder names
confounders <- c(grep("^site_", colnames(confounds_dummy), value=TRUE),
                 'sex', 'age', 'date_diff',
                 'head_size', 'head_location1', 'head_location2','head_location3','head_location4',
                 'qc_discrepancy', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5')
sites <- grep("^site_", colnames(confounds_dummy), value=TRUE)
head_locations <- c('head_location1', 'head_location2','head_location3','head_location4')
genes <- c('gene1', 'gene2', 'gene3', 'gene4', 'gene5')

## Select confounder sets
#only age
ukb_z1 <- confounds_dummy[, .SD, .SDcols = c('id', 'age')]
data.table::fwrite(ukb_z1, file = paste0(path_to_save_confounders, "ukb_z1_age.csv"))
#age, head size
ukb_z2 <- confounds_dummy[, .SD, .SDcols = c('id', 'age', 'head_size')]
data.table::fwrite(ukb_z2, file = paste0(path_to_save_confounders, "ukb_z2_agesex.csv"))
#age, sex, head size, site
ukb_z4 <- confounds_dummy[, .SD, .SDcols = c('id', 'age', 'head_size', sites, 'sex')]
data.table::fwrite(ukb_z4, file = paste0(path_to_save_confounders, "ukb_z4_agesexsitesize.csv"))
#age, sex, head size, site, date, qc-discrepancy
ukb_z6 <- confounds_dummy[, .SD, .SDcols = c('id', 'age', 'head_size', sites, 'sex', 'date_diff', 'qc_discrepancy')]
data.table::fwrite(ukb_z6, file = paste0(path_to_save_confounders, "ukb_z6_agesexsitesizedateqc.csv"))
#age, sex, head size, site, date, qc-discrepancy, 4xhead location
ukb_z10 <- confounds_dummy[, .SD, .SDcols = c('id', 'age', 'head_size', sites, 'sex', 'date_diff', 'qc_discrepancy', head_locations)]
data.table::fwrite(ukb_z10, file = paste0(path_to_save_confounders, "ukb_z10_agesexsitesizedateqclocation.csv"))
#age, sex, head size, site, date, qc-discrepancy, 4xhead location, 5xgenes
ukb_z15 <- confounds_dummy[, .SD, .SDcols = c('id', 'age', 'head_size', sites, 'sex', 'date_diff', 'qc_discrepancy', head_locations, genes)]
data.table::fwrite(ukb_z15, file = paste0(path_to_save_confounders, "ukb_z15_agesexsitesizedateqcgenes.csv"))
