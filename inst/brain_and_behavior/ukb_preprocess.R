##### Preprocessing of UKB data for conditional associations between brain structures and behavioral traits
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(paletteer)
library('DNCIT')
library(fastDummies)


####### Real-world application
# paths to data repository
path_to_ukb <- ''
path_to_ukb_data <- paste0(path_to_ukb, 'ukb49727.csv')
path_to_fastsurfer_ids <- paste0(path_to_ukb, 'IDPs/Ids_IDPs.csv')
path_to_fractional_anisotropy_ids <- paste0(path_to_ukb, 'Real-world_application/ids_FA_measures.csv')
path_to_freesurfer_dk_atlas_ids <- paste0(path_to_ukb, 'Real-world_application/ids_freesurfer_Desikan-Killiany_atlas.csv')
path_to_freesurfer_aseg_ids <- paste0(path_to_ukb, 'Real-world_application/ids_freesurfer_ASEG.csv')
path_to_save_ids_brain_avinun <- paste0(path_to_ukb, 'ids/ids_brain_avinun.csv')
path_to_save_ids_personality <- paste0(path_to_ukb, 'ids/ids_personality.csv')
path_to_save_ids_confounder_avinun <- paste0(path_to_ukb, 'ids/ids_confounder_avinun.csv')
path_to_save_ids_personality_avinun <- paste0(path_to_ukb, 'ids/ids_personality_avinun.csv')
# paths to save results
path_to_save_preprocessed_data <- paste0(path_to_ukb, 'ukb_free_fast_behavior_healthy.csv')

##color palettes for plotting
palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")

### Subset available UKB to obtain available col names
ukb_whole_columns <- data.table::fread(file=path_to_ukb_data, header=TRUE, nrows=10)

##### behavioural traits and MRIs
### 1) Load data based on ids
## general ids
id_eid = 'eid'

### ids to remove subjects
# cancer (sum over all visits)
ids_cancer <- colnames(ukb_whole_columns)[grepl('2453', colnames(ukb_whole_columns), fixed = TRUE)]
## ICD10 (diagnosis)
ids_ICD10 <- colnames(ukb_whole_columns)[grepl('41202', colnames(ukb_whole_columns), fixed = TRUE)]
## Use of psychotropic, glucocorticoid, or hypolipidemic medication based on medication
ids_medication <- colnames(ukb_whole_columns)[grepl('20003', colnames(ukb_whole_columns), fixed = TRUE)]
#all
ids_remove <- c(ids_cancer, ids_ICD10, ids_medication)

## ids brain morphometry
# Freesurfer and Fastsurfer
ids_IDPs <- data.table::fread(path_to_fastsurfer_ids, select='Field ID')
ids_IDPs_full_str <- ids_IDPs %>%
  dplyr::mutate_all(list(~ stringr::str_c(., '-2.0'))) %>%
  dplyr::pull(`Field ID`)
# fractional anisotropy
ids_FA_lr <- stats::na.omit(data.table::fread(file=path_to_fractional_anisotropy_ids,header=TRUE))
ids_FA <- rep(0, nrow(ids_FA_lr))
for (i in 1:nrow(ids_FA_lr)){
  ids_FA[i] <- colnames(ukb_whole_columns)[grepl(paste(ids_FA_lr[i,1],'-2.0', sep=""), colnames(ukb_whole_columns), fixed = TRUE)]
}
ids_brain_morphometry <- c(ids_IDPs_full_str,ids_FA)

## ids big five behavioural traits
# Neuroticism score (derived in Smith et al. 2013)
id_neuroticism <- colnames(ukb_whole_columns)[grepl('20127', colnames(ukb_whole_columns), fixed = TRUE)]
# Sociability
id_friends_visits <- colnames(ukb_whole_columns)[grepl('1031-2.0', colnames(ukb_whole_columns), fixed = TRUE)]
id_guilty_feelings <- colnames(ukb_whole_columns)[grepl('\\<2030-2.0\\>', colnames(ukb_whole_columns))]
id_tired <- colnames(ukb_whole_columns)[grepl('\\<2080-2.0\\>', colnames(ukb_whole_columns))]
id_leisure <- colnames(ukb_whole_columns)[grepl('\\<6160-2.0\\>', colnames(ukb_whole_columns))]
ids_social <- c(id_friends_visits, id_guilty_feelings, id_tired, id_leisure)
# warmth
id_confide <- colnames(ukb_whole_columns)[grepl('\\<2110-2.0\\>', colnames(ukb_whole_columns))]
id_irritable <- colnames(ukb_whole_columns)[grepl('\\<1940-2.0\\>', colnames(ukb_whole_columns))]
id_moody <- colnames(ukb_whole_columns)[grepl('\\<1920-2.0\\>', colnames(ukb_whole_columns))]
id_tense <- colnames(ukb_whole_columns)[grepl('\\<1990-2.0\\>', colnames(ukb_whole_columns))]
id_nervous <- colnames(ukb_whole_columns)[grepl('\\<1970-2.0\\>', colnames(ukb_whole_columns))]
ids_warmth <- c(id_confide, id_irritable, id_moody, id_tense, id_nervous)
# diligence
id_no_interest <- colnames(ukb_whole_columns)[grepl('\\<2060-2.0\\>', colnames(ukb_whole_columns))]
id_fed_up <- colnames(ukb_whole_columns)[grepl('\\<1960-2.0\\>', colnames(ukb_whole_columns))]
id_risk <- colnames(ukb_whole_columns)[grepl('\\<2040-2.0\\>', colnames(ukb_whole_columns))]
id_worry_long <- colnames(ukb_whole_columns)[grepl('\\<2000-2.0\\>', colnames(ukb_whole_columns))]
ids_diligence <- c(id_no_interest, id_fed_up, id_risk, id_worry_long)
# curiosity
id_lonely <- colnames(ukb_whole_columns)[grepl('\\<2020-2.0\\>', colnames(ukb_whole_columns))]
id_nerves <- colnames(ukb_whole_columns)[grepl('\\<2010-2.0\\>', colnames(ukb_whole_columns))]
id_tense_last_two_weeks <- colnames(ukb_whole_columns)[grepl('\\<2070-2.0\\>', colnames(ukb_whole_columns))]
ids_curiosity <- c(id_lonely, id_nerves, id_tense_last_two_weeks, id_risk)
#nervousness
id_easy_hurt <- colnames(ukb_whole_columns)[grepl('\\<1950-2.0\\>', colnames(ukb_whole_columns))]
ids_nervousness <- c(id_tense, id_irritable, id_no_interest, id_moody, id_easy_hurt)
# all
ids_personality <- c(id_neuroticism,
                     id_friends_visits,id_guilty_feelings, id_tired, id_leisure,
                     id_confide, id_irritable, id_moody, id_tense, id_nervous,
                     id_no_interest, id_fed_up, id_risk, id_worry_long,
                     id_lonely, id_nerves, id_tense_last_two_weeks,
                     id_easy_hurt)
#data.table::fwrite(as.data.frame(ids_personality), file = path_to_save_ids_personality)

## ids confounders
# age, sex, assessment center, qc discrepency, aquisation date, head size
id_sex <- '31-0.0'
id_age_assessment_center <- '21003-2.0'
id_assessment_center <- '54-0.0'
id_home_location_urban <- '20118-0.0'
id_qc_discrepancy <- "25732-2.0"
id_date <- "53-2.0"
id_head_size <- "25000-2.0"
# location head in scanner
ids_head_location <- c("25756-2.0", "25757-2.0", "25758-2.0", "25759-2.0")
# T2_flair used for T1

# genetics
ids_genetic_pcs <- colnames(ukb_whole_columns)[grepl('22009', colnames(ukb_whole_columns), fixed = TRUE)][1:5]
# head motion (for fractional anisotropy)
#id_head_motion <- colnames(ukb_whole_columns)[grepl('4244-1.9', colnames(ukb_whole_columns), fixed = TRUE)]
#ukb_data <- data.table::fread(file=paste0(path_to_ukb, 'ukb49727.csv'),select = id_head_motion,  header=TRUE, nrows=100)
ids_confounders <- c(id_sex, id_age_assessment_center,
                     id_assessment_center, id_date,
                     id_head_size, ids_head_location, id_qc_discrepancy,
                     ids_genetic_pcs)
#data.table::fwrite(as.data.frame(ids_confounders), file = path_to_save_ids_confounder_avinun)

#### Full data
ids_ukb_mri_personality <- c(id_eid,
                             ids_confounders,
                             ids_remove,
                             ids_brain_morphometry,
                             ids_personality)
#data.table::fwrite(as.data.frame(ids_ukb_mri_personality), file = path_to_save_ids_personality_avinun)
# Read data depending on selected ids
ukb_data <- data.table::fread(file=path_to_ukb_data,select = ids_ukb_mri_personality,  header=TRUE)
no_MRIs_fastsurfer1 <- !apply(dplyr::select(ukb_data, dplyr::all_of(ids_IDPs_full_str[1:2])),1, anyNA)
ukb_data_mri_1 <- ukb_data[no_MRIs_fastsurfer1,]
no_MRIs_fastsurfer <- !apply(dplyr::select(ukb_data_mri_1, dplyr::all_of(ids_IDPs_full_str)),1, anyNA)
ukb_data_mri <- ukb_data_mri_1[no_MRIs_fastsurfer,]



### 2) Data preprocessingRemove subjects with:
## 1) medical diagnoses of cancer, stroke, diabetes requiring insulin treatment, chronic kidney or liver disease, or lifetime history of psychotic symptoms
#cancer (sum over all visits)
ukb_data_cancer <- dplyr::mutate_all(dplyr::select(ukb_data_mri, dplyr::all_of(ids_cancer)), as.character)
patients_w_cancer <- apply(ukb_data_cancer, 1, function(row) {
  any(sapply('1', function(substring) {
    any(grepl(substring, row))
  }))
})
## ICD10
ids_ICD10 <- colnames(ukb_whole_columns)[grepl('41202', colnames(ukb_whole_columns), fixed = TRUE)]
ukb_data_ICD10 <- dplyr::select(ukb_data_mri, dplyr::all_of(ids_ICD10))
# stroke
ICD10_stroke <- 'I64'
patients_w_stroke <- apply(ukb_data_ICD10, 1, function(row) {
  any(sapply(row, function(x) grepl(ICD10_stroke, as.character(x))))
})
#insulin (substring check)
ICD10_diabetes_w_insulin <- 'E10'
patients_w_diabetes <- apply(ukb_data_ICD10, 1, function(row) {
  any(sapply(row, function(x) grepl(ICD10_diabetes_w_insulin, as.character(x))))
})
#kidney disease
ICD10_kidney_failure <- c('E883', 'I120', 'I129', 'I1310', 'I1311', 'R34', 'T795XXA', 'N17', 'N18', 'N19')
patients_w_kidney_disease <- apply(ukb_data_ICD10, 1, function(row) {
  any(sapply(ICD10_kidney_failure, function(substring) {
    any(grepl(substring, row))
  }))
})
# liver disease
ICD10_liver_failure <- 'K7'
patients_w_liver_disease <- apply(ukb_data_ICD10, 1, function(row) {
  any(sapply(row, function(x) grepl(ICD10_liver_failure, as.character(x))))
})
#no personal history of mental and behavioural disorders
ICD10_history_mental <- 'Z865'
patients_w_history_mental <- apply(ukb_data_ICD10, 1, function(row) {
  any(sapply(row, function(x) grepl(ICD10_history_mental, as.character(x))))
})
## 2) use of psychotropic, glucocorticoid, or hypolipidemic (missing currently) medication
ukb_data_meds <- dplyr::select(ukb_data_mri, dplyr::all_of(ids_medication))
# psychotropic
psychotropic_meds_drug_codes <- data.table::fread(file=paste0(path_to_ukb, 'Real-world_application/psychotropic_medication_list_codes.csv'),select='Drug code',header=TRUE)
psychotropic_meds_drug_codes <- unlist(dplyr::mutate_all(psychotropic_meds_drug_codes, as.character))
patients_w_psychotropic_meds <- apply(ukb_data_meds, 1, function(row) {
  any(sapply(psychotropic_meds_drug_codes, function(substring) {
    any(grepl(substring, row))
  }))
})
# glucocorticoid
glucocorticoid_meds_drug_codes1 <- data.table::fread(file=paste0(path_to_ukb, 'Real-world_application/glucocorticoids_meds_list_codes.csv'),select='Code',header=TRUE)
glucocorticoid_meds_drug_codes2 <- data.table::fread(file=paste0(path_to_ukb, 'Real-world_application/glucocorticoids__meds_list_codes2.csv'),select='Code',header=TRUE)
glucocorticoid_meds_drug_codes <- unlist(dplyr::mutate_all(rbind(glucocorticoid_meds_drug_codes1, glucocorticoid_meds_drug_codes2), as.character))
patients_w_glucocorticoid_meds <- apply(ukb_data_meds, 1, function(row) {
  any(sapply(glucocorticoid_meds_drug_codes, function(substring) {
    any(grepl(substring, row))
  }))
})
## Additionally 4) personality and psychiatric disorders
ICD10_personality_psychiatric_disorders <- c('F07', 'F09', 'F1', 'F6', 'F2', 'F3', 'F4', 'F50')
patients_w_pers_psy_disorder <- apply(ukb_data_ICD10, 1, function(row) {
  any(sapply(ICD10_personality_psychiatric_disorders, function(substring) {
    any(grepl(substring, row))
  }))
})

## Remove if any of the above is TRUE for individual
remove_vectors <- list(
  patients_w_cancer, patients_w_stroke, patients_w_diabetes, patients_w_kidney_disease,patients_w_liver_disease,
  patients_w_history_mental, patients_w_psychotropic_meds, patients_w_glucocorticoid_meds, patients_w_pers_psy_disorder
)
remove_binary <- apply(simplify2array(remove_vectors), 1, any)
subset_ukb_data_mri <- as.data.frame(ukb_data_mri[remove_binary,])


####Replication Avinun Analysis UK Biobank
## brain structure Freesurfer measurement preparation
ids_freesurfer_atlas <- data.table::fread(file=path_to_freesurfer_dk_atlas_ids,header=TRUE)
ids_freesurfer_aseg <- stats::na.omit(data.table::fread(file=path_to_freesurfer_aseg_ids,header=TRUE))
ids_FA_all <- stats::na.omit(data.table::fread(file=path_to_fractional_anisotropy_ids,header=TRUE))
ids_FA_lr <- ids_FA_all[grepl('left', ids_FA_all$Description) | grepl('right', ids_FA_all$Description) ,]
ids_surface_area_lr <- ids_freesurfer_atlas[1:68,]
ids_cortical_thick_lr <- ids_freesurfer_atlas[69:136,]
ids_gray_matter_vols_lr <- ids_freesurfer_aseg[c(47:50,61:64,67:68,72:73,79:82, 86:87),]
ids_lr <- rbind(ids_cortical_thick_lr,ids_surface_area_lr, ids_gray_matter_vols_lr, ids_FA_lr)
## Average measures of SA, CT and extra
for(i in 1:(nrow(ids_lr)/2)){
  print(i)
  ids_select <- c(colnames(subset_ukb_data_mri)[grepl(ids_lr[2*i-1,1], colnames(subset_ukb_data_mri))],
                  colnames(subset_ukb_data_mri)[grepl(ids_lr[2*i,1], colnames(subset_ukb_data_mri))])
  measurements <- subset_ukb_data_mri[, ids_select]
  #average and scale brain measurements
  avg <- scale(rowMeans(measurements))
  new_id <- paste(as.character(unlist(ids_lr[2*i-1,1])))
  print(new_id)
  subset_ukb_data_mri[[new_id]] <- avg
  #print(subset_ukb_data_mri[1:5,new_id])
}
## whole brain measurements
ids_whole_brain_aseg <- ids_freesurfer_aseg[c(71,83, 88)]
ids_FA_whole <- ids_FA_all[ids_FA_all$ID %notin% ids_FA_lr$ID ,]
ids_whole_brain_aseg <- rbind(ids_whole_brain_aseg, ids_FA_whole)
ids_whole_brain_aseg_vec <- rep(1, nrow(ids_whole_brain_aseg))
for(i in 1:nrow(ids_whole_brain_aseg)){
  ids_whole_brain_aseg_vec[i] <- colnames(subset_ukb_data_mri)[grepl(ids_whole_brain_aseg[i,1], colnames(subset_ukb_data_mri))]
  #scale
  subset_ukb_data_mri[, ids_whole_brain_aseg_vec[i]] <- scale(subset_ukb_data_mri[ids_whole_brain_aseg_vec[i]])
}

## save subset
#data.table::fwrite(subset_ukb_data_mri, file = path_to_save_preprocessed_data)
ids_brain_structure <- c(as.character(unlist(ids_lr[seq(from = 1, to = nrow(ids_lr), by = 2),1])),
                         ids_whole_brain_aseg_vec)
ids_brain <- rbind(ids_lr[seq(from = 1, to = nrow(ids_lr), by = 2),], ids_whole_brain_aseg)
brain_structures <- paste0("brain_structure_", 1:length(ids_brain_structure))
## save ids brain structures in Avinun 2020
#data.table::fwrite(as.data.frame(ids_brain_structure), file = path_to_save_ids_brain_avinun)
