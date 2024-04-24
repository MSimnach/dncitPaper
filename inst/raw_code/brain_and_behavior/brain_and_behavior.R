##### Analysis of conditional associations between brain structures and behavioral traits in the UKB cohort
library(data.table)
#library(tidyverse)
library(dplyr)
#library(lattice)
#library(ggcorrplot)
#setwd('C:\\Users\\Marco\\seadrive_root\\Marco Si\\Meine Bibliotheken\\Meine Bibliothek\\CITs\\P1 Deep CITs\\UK biobank\\Data\\IDPs\\vols')

#CITs
#devtools::install('C:/Users/Marco/seadrive_root/Marco Si/Meine Bibliotheken/Meine Bibliothek/Coding/DNCIT')
library('DNCIT')
#library(lmtest)
#library(sandwich)
#library(caret)
#library(fastDummies)


####### Real-world application
### Subset available UKB to obtain available col names
ukb_whole_columns <- data.table::fread(file='M:/CITs/Application/UKB_data/ukb49727.csv', header=TRUE, nrows=10)

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
ids_IDPs <- fread('C:\\Users\\Marco\\seadrive_root\\Marco Si\\Meine Bibliotheken\\Meine Bibliothek\\CITs\\P1 Deep CITs\\UK biobank\\Data/IDPs/Ids_IDPs.csv', select='Field ID')
ids_IDPs_full_str <- ids_IDPs %>%
  mutate_all(list(~ str_c(., '-2.0'))) %>%
  pull(`Field ID`)
# fractional anisotropy
ids_FA_lr <- na.omit(fread(file='M:/CITs/Application/UKB_data/Real-world_application/ids_FA_measures.csv',header=TRUE))
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
#ukb_data <- fread(file='M:/CITs/Application/UKB_data/ukb49727.csv',select = id_head_motion,  header=TRUE, nrows=100)
ids_confounders <- c(id_sex, id_age_assessment_center,
                     id_assessment_center, id_date,
                     id_head_size, ids_head_location, id_qc_discrepancy,
                     ids_genetic_pcs)

#### Full data
ids_ukb_mri_personality <- c(id_eid,
                             ids_confounders,
                             ids_remove,
                             ids_brain_morphometry,
                             ids_personality)
# Read data depending on selected ids
ukb_data <- fread(file='M:/CITs/Application/UKB_data/ukb49727.csv',select = ids_ukb_mri_personality,  header=TRUE)
no_MRIs_fastsurfer1 <- !apply(select(ukb_data, all_of(ids_IDPs_full_str[1:2])),1, anyNA)
ukb_data_mri_1 <- ukb_data[no_MRIs_fastsurfer1,]
no_MRIs_fastsurfer <- !apply(select(ukb_data_mri_1, all_of(ids_IDPs_full_str)),1, anyNA)
ukb_data_mri <- ukb_data_mri_1[no_MRIs_fastsurfer,]


### 2) Data preprocessingRemove subjects with:
## 1) medical diagnoses of cancer, stroke, diabetes requiring insulin treatment, chronic kidney or liver disease, or lifetime history of psychotic symptoms
#cancer (sum over all visits)
ukb_data_cancer <- mutate_all(select(ukb_data_mri, all_of(ids_cancer)), as.character)
patients_w_cancer <- apply(ukb_data_cancer, 1, function(row) {
  any(sapply('1', function(substring) {
    any(grepl(substring, row))
  }))
})
## ICD10
ids_ICD10 <- colnames(ukb_whole_columns)[grepl('41202', colnames(ukb_whole_columns), fixed = TRUE)]
ukb_data_ICD10 <- select(ukb_data_mri, all_of(ids_ICD10))
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
## 2) use of psychotropic, glucocorticoid, or hypolipidemic medication
ukb_data_meds <- select(ukb_data_mri, all_of(ids_medication))
# psychotropic
psychotropic_meds_drug_codes <- fread(file='M:/CITs/Application/UKB_data/Real-world_application/psychotropic_medication_list_codes.csv',select='Drug code',header=TRUE)
psychotropic_meds_drug_codes <- unlist(mutate_all(psychotropic_meds_drug_codes, as.character))
patients_w_psychotropic_meds <- apply(ukb_data_meds, 1, function(row) {
  any(sapply(psychotropic_meds_drug_codes, function(substring) {
    any(grepl(substring, row))
  }))
})
# glucocorticoid
glucocorticoid_meds_drug_codes1 <- fread(file='M:/CITs/Application/UKB_data/Real-world_application/glucocorticoids_meds_list_codes.csv',select='Code',header=TRUE)
glucocorticoid_meds_drug_codes2 <- fread(file='M:/CITs/Application/UKB_data/Real-world_application/glucocorticoids__meds_list_codes2.csv',select='Code',header=TRUE)
glucocorticoid_meds_drug_codes <- unlist(mutate_all(rbind(glucocorticoid_meds_drug_codes1, glucocorticoid_meds_drug_codes2), as.character))
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
ids_freesurfer_atlas <- fread(file='M:/CITs/Application/UKB_data/Real-world_application/ids_freesurfer_Desikan-Killiany_atlas.csv',header=TRUE)
ids_freesurfer_aseg <- na.omit(fread(file='M:/CITs/Application/UKB_data/Real-world_application/ids_freesurfer_ASEG.csv',header=TRUE))
ids_FA_lr <- na.omit(fread(file='M:/CITs/Application/UKB_data/Real-world_application/ids_FA_measures.csv',header=TRUE))
ids_cortical_thick_lr <- ids_freesurfer_atlas[1:68,]
ids_surface_area_lr <- ids_freesurfer_atlas[69:136,]
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
  print(subset_ukb_data_mri[1:5,new_id])
}
## whole brain measurements
ids_whole_brain_aseg <- ids_freesurfer_aseg[c(71,83, 88)]
ids_whole_brain_aseg_vec <- rep(1, nrow(ids_whole_brain_aseg))
for(i in 1:nrow(ids_whole_brain_aseg)){
  ids_whole_brain_aseg_vec[i] <- colnames(subset_ukb_data_mri)[grepl(ids_whole_brain_aseg[i,1], colnames(subset_ukb_data_mri))]
  #scale
  subset_ukb_data_mri[, ids_whole_brain_aseg_vec[i]] <- scale(subset_ukb_data_mri[, ids_whole_brain_aseg_vec[i]])
}
ids_brain_structure <- c(as.character(unlist(ids_lr[seq(from = 1, to = nrow(ids_lr), by = 2),1])),
                         ids_whole_brain_aseg_vec)
brain_structures <- paste0("brain_structure_", 1:107)


#### Analysis
##Variables used and data subset
ids_ukb_mri_neuroticism_analyis <- c(ids_confounders, ids_brain_structure, unique(ids_personality))
mri_neuroticism <- na.omit(subset_ukb_data_mri[,ids_ukb_mri_neuroticism_analyis])

###Outputs
# Neuroticism
mri_neuroticism$neuroticism <- scale(mri_neuroticism[, id_neuroticism])
# Sociability score
mri_neuroticism$sociability <- scale(rowSums(mri_neuroticism[, ids_social]))
# Warmth score
mri_neuroticism$warmth <- scale(rowSums(mri_neuroticism[, ids_warmth]))
# Diligence score
mri_neuroticism$diligence <- scale(rowSums(mri_neuroticism[, ids_diligence]))
# Curiosity score
mri_neuroticism$curiosity <- scale(rowSums(mri_neuroticism[, ids_curiosity]))
# Nervousness score
mri_neuroticism$nervousness <- scale(rowSums(mri_neuroticism[, ids_nervousness]))

## remove all original columns used to obtain the scores
mri_neuroticism <- mri_neuroticism[, !names(mri_neuroticism) %in% unique(ids_personality)]
personality_traits <- c('neuroticism', 'sociability', 'warmth', 'diligence', 'curiosity', 'nervousness')
colnames(mri_neuroticism) <- c('sex', 'age', 'site', 'date',
                               'head_size', 'head_location1', 'head_location2','head_location3','head_location4',
                               'qc_discrepancy', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5',
                               brain_structures, personality_traits)
## dummy confounder site
dummy_df <- dummy_cols(mri_neuroticism, select_columns = 'site')
mri_neuroticism_dummy <- dummy_df[, !names(dummy_df) %in% c('site', 'site_10003')]
## date difference to first date
mri_neuroticism_dummy$date_diff <- mri_neuroticism_dummy[,'date'] - min(mri_neuroticism_dummy[,'date'])
mri_neuroticism_dummy <- mri_neuroticism_dummy[, !names(mri_neuroticism_dummy) %in% c('date')]

# confounder names
confounders <- c(grep("^site_", colnames(mri_neuroticism_dummy), value=TRUE),
                 'sex', 'age', 'date_diff',
                 'head_size', 'head_location1', 'head_location2','head_location3','head_location4',
                 'qc_discrepancy', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5')

### Coefficient tests for each brain measurement and trait in linear model with robust standard errors and adj. p-values
lm_models <- list()
#pvalues_lrtest_neuroticism <- list()
pvalues_ttest <- list()
# Loop over each brain structure
for (brain_structure in brain_structures) {
  print(brain_structure)
  #regress out head location in scanner
  for (trait in personality_traits) {
    #df for regression
    df_brain_structure <- mri_neuroticism_dummy[,c(confounders, brain_structure, trait)]

    # Fit linear model for the current brain structure
    lm_formula <- as.formula(paste(brain_structure, "~ . + I(age^2)+age:sex+I(date_diff^2) -", brain_structure))
    lm_model <- lm(lm_formula, data = df_brain_structure)

    # Store the linear model in the list
    lm_models[[paste0(brain_structure, trait)]] <- lm_model
    # Perform coefficient test for trait
    #pvalue_lrtest <- lrtest(lm_model, trait)
    pvalue_ttest <- coeftest(lm_model, vcov = vcovHC(lm_model, type='HC0'), parm=trait)[trait,4]

    # Store the test result in the list
    pvalues_ttest[[paste0(brain_structure, trait)]] <- pvalue_ttest
    #pvalues_lrtest_neuroticism[[paste0(brain_structure, trait)]] <- pvalue_lrtest$`Pr(>Chisq)`[2]
  }
}
## p-values
p_val_df <- do.call(rbind, lapply(pvalues_ttest, as.data.frame))
p_val_df$p_adj <- p.adjust(p_val_df[,1], method='BH')
colnames(p_val_df) <- c('p_unadj', 'p_adj')
p_val_df[order(p_val_df$p_unadj),]

#### joint significance of brain structures
## Wald test for each trait separately
lm_models_joint <- list()
p_joint_wald <- list()
for (trait in personality_traits){
  lm_formula <- as.formula(paste(trait, "~ . + I(age^2)+age:sex+I(date_diff^2) -", trait))
  lm_model <- lm(lm_formula, mri_neuroticism_dummy[, !names(mri_neuroticism_dummy) %in% personality_traits[!(personality_traits %in% trait)]])
  lm_models_joint[[trait]] <- lm_model
  p_wald <- waldtest(lm_model, brain_structures)
  p_joint_wald[[trait]] <- p_wald$`Pr(>F)`[2]
}
p_joint_df <- do.call(rbind, lapply(p_joint_wald, as.data.frame))
p_joint_df$p_adj <- rep(1, length(personality_traits))
colnames(p_joint_df) <- c('p_unadj', 'p_adj')
p_total_wald <- rbind(p_val_df, p_joint_df)
p_total_wald$p_adj <- p.adjust(p_total_wald[,1], method='BH')
p_total_wald[order(p_total_wald$p_unadj),]

### RCoT
rcots_p_joint_traits <- list()
# for each trait
for (trait in personality_traits){
  print(trait)
  X <- as.matrix(mri_neuroticism_dummy[, brain_structures])
  Y <- as.matrix(mri_neuroticism_dummy[, trait])
  Z <- sapply(as.data.frame(mri_neuroticism_dummy[, confounders]), as.numeric)
  rcots_p_joint_trait <- rep(1,20)
  for(seed in 1:20){
    rcot_p <- RCoT(X,Y,Z, seed = seed)$p
    rcots_p_joint_trait[seed] <- rcot_p
  }
  rcots_p_joint_traits[[trait]] <- mean(rcots_p_joint_trait)
}
# for all traits together
X <- as.matrix(mri_neuroticism_dummy[, brain_structures])
Y <- as.matrix(mri_neuroticism_dummy[, personality_traits])
Z <- sapply(as.data.frame(mri_neuroticism_dummy[, confounders]), as.numeric)
rcots_p_joint_all_traits <- rep(1,20)
for(seed in 1:20){
  rcot_p <- RCoT(X,Y,Z, seed = seed)$p
  rcots_p_joint_all_traits[seed] <- rcot_p
}
rcots_p_joint_traits[['all_traits']] <- mean(rcots_p_joint_all_traits)



#### Confounders and DNCITs
### 1) check if nonlinear relationships with confounders are satisfactory addressed
confounders_short <- c('age', 'sex', 'genes', 'head_size', 'head_position', 'date_diff', 'site')
brain_structures_lin_unconfounded <- list()
# for each confounder compute linearly unconfounded brain structure measurements
for(confounder in confounders_short){
  count <- 1
  brain_structures_lin_unconfounded[[confounder]] <- matrix(nrow = nrow(mri_neuroticism_dummy), ncol=length(brain_structures))
  for(brain_structure in brain_structures){
    if (confounder == 'head_position'){
      lm_formula <- as.formula(paste(brain_structure,'~head_location1+head_location2+head_location3+head_location4'))
    }else if (confounder == 'genes'){
      lm_formula <- as.formula(paste(brain_structure,'~gene1+gene2+gene3+gene4+gene5'))
    }else if (confounder == 'site'){
      lm_formula <- as.formula(paste(brain_structure,'~', paste(grep("^site_", colnames(mri_neuroticism_dummy), value=TRUE), collapse = "+")))
    }else{
      lm_formula <- as.formula(paste(brain_structure,'~', confounder))
    }

    lm_model <- lm(lm_formula, mri_neuroticism_dummy)
    brain_structures_lin_unconfounded[[confounder]][, count] <- lm_model$residuals
    count <- count +1
  }
  colnames(brain_structures_lin_unconfounded[[confounder]]) <- brain_structures
}
# CIT unconfounded residuals \indep confounder | other_confounders
mri_lin_unconf <- list()
rcots_p_mean <- list()
for (confounder in confounders_short){
  print(confounder)
  mri_lin_unconf[[confounder]] <- cbind(mri_neuroticism_dummy[, !names(mri_neuroticism_dummy) %in% brain_structures], brain_structures_lin_unconfounded[[confounder]])
  X <- as.matrix(mri_lin_unconf[[confounder]][, brain_structures])
  if (confounder == 'head_position'){
    head_pos <- c('head_location1','head_location2','head_location3','head_location4')
    confound_var <- as.matrix(mri_lin_unconf[[confounder]][, head_pos])
    Z <- sapply(as.data.frame(mri_lin_unconf[[confounder]][, confounders %notin% head_pos]), as.numeric)
  }else if (confounder == 'genes'){
    genes_ <- c('gene1','gene2','gene3','gene4','gene5')
    confound_var <- as.matrix(mri_lin_unconf[[confounder]][, genes_])
    Z <- sapply(as.data.frame(mri_lin_unconf[[confounder]][, confounders %notin% genes_]), as.numeric)
  }else if (confounder == 'site'){
    sites_ <- grep("^site_", colnames(mri_neuroticism_dummy), value=TRUE)
    confound_var <- as.matrix(mri_lin_unconf[[confounder]][, sites_])
    Z <- sapply(as.data.frame(mri_lin_unconf[[confounder]][, confounders %notin% sites_]), as.numeric)
  }else{
    confound_var <- as.matrix(mri_lin_unconf[[confounder]][, confounder])
    Z <- sapply(as.data.frame(mri_lin_unconf[[confounder]][, confounders %notin% confounder]), as.numeric)
  }
  rcots_p <- rep(1, 20)
  for(seed in 1:20){
    rcots_p[seed] <- RCoT(X,Y,Z, seed = seed)$p
  }
  rcots_p_mean[[confounder]] <- mean(rcots_p)
}
# still evidence for nonlinear confounding for age, sex, genes, head size, head position and date diff




### 2) Additional confounders
n
X <- as.matrix(mri_neuroticism_dummy[, brain_structures])
confounder <- as.matrix(mri_neuroticism_dummy[, personality_traits])
Z <- sapply(as.data.frame(mri_neuroticism_dummy[, names(mri_neuroticism_dummy) %in% confounders[confounders!='age']]), as.numeric)
age <- mri_neuroticism_dummy[,'age']
