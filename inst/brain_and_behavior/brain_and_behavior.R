##### Analysis of conditional associations between brain structures and behavioral traits in the UKB cohort
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
path_to_ukb_data <- paste0(path_to_ukb,'ukb49727.csv')
path_to_fastsurfer_ids <- paste0(path_to_ukb,'ids/ids_IDPs.csv')
path_to_fractional_anisotropy_ids <- paste0(path_to_ukb,'ids/ids_FA_measures.csv')
path_to_freesurfer_dk_atlas_ids <- paste0(path_to_ukb,'ids/ids_freesurfer_Desikan-Killiany_atlas.csv')
path_to_freesurfer_aseg_ids <- paste0(path_to_ukb,'ids/ids_freesurfer_ASEG.csv')
path_to_save_ids_brain_avinun <- paste0(path_to_ukb,'ids/ids_brain_avinun.csv')
path_to_save_ids_confounder_avinun <- paste0(path_to_ukb,'ids/ids_confounder_avinun.csv')
path_to_save_ids_personality_avinun <- paste0(path_to_ukb,'ids/ids_personality_avinun.csv')
path_to_save_ids_personality <- paste0(path_to_ukb,'ids/ids_personality.csv')
# paths to save results
path_to_save_preprocessed_data <- paste0(path_to_ukb,'ukb_free_fast_behavior_healthy.csv')
path_to_save_ukb_avinun <- paste0(path_to_ukb,'ukb_avinun.csv')
path_to_save_ukb_fast_behavior <- paste0(path_to_ukb,'ukb_fast_behavior.csv')
path_to_pval_brain_trait_each <- paste0(path_to_ukb,'Real-world_application/p_val_structure_trait.csv')
path_to_pval_brain_trait_each_w_joint <- paste0(path_to_ukb,'Real-world_application/p_val_structure_trait_w_joint_structures.csv')
path_to_pval_brain_trait_joint_rcot <- paste0(path_to_ukb,'Real-world_application/p_val_structures_trait.csv')
path_to_pval_brain_trait_confounder_control <- paste0(path_to_ukb,'Real-world_application/p_val_structures_trait_confounder_control.csv')
path_to_save_plots <- ''

##color palettes for plotting
#palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")
palet_discrete <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")

### Subset available UKB to obtain available col names
ukb_whole_columns <- data.table::fread(file=path_to_ukb_data, header=TRUE, nrows=10)


#### Analysis
## load preprocessed data and ids
ids_brain_structure <- data.table::fread(file=path_to_save_ids_brain_avinun,header=TRUE)$ids_brain_structure
brain_structures <- paste0("brain_structure_", 1:length(ids_brain_structure))
ids_confounders <- data.table::fread(file=path_to_save_ids_confounder_avinun,header=TRUE)$ids_confounder
ids_personality <- data.table::fread(file=path_to_save_ids_personality,header=TRUE)$ids_personality
subset_ukb_data_mri <- as.data.frame(data.table::fread(file=path_to_save_preprocessed_data,header=TRUE))
##Variables used and data subset
ids_ukb_mri_neuroticism_analyis <- c(ids_confounders, ids_brain_structure, ids_personality)
mri_neuroticism <- stats::na.omit(subset_ukb_data_mri[,ids_ukb_mri_neuroticism_analyis])

##Standardize Y
#ids for Y
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
dummy_df <- fastDummies::dummy_cols(mri_neuroticism, select_columns = 'site')
mri_neuroticism_dummy <- dummy_df[, !names(dummy_df) %in% c('site', 'site_10003')]
## date difference to first date
mri_neuroticism_dummy$date_diff <- mri_neuroticism_dummy[,'date'] - min(mri_neuroticism_dummy[,'date'])
mri_neuroticism_dummy <- mri_neuroticism_dummy[, !names(mri_neuroticism_dummy) %in% c('date')]
#data.table::fwrite(mri_neuroticism_dummy, file = path_to_save_ukb_avinun)

# confounder names
confounders <- c(grep("^site_", colnames(mri_neuroticism_dummy), value=TRUE),
                 'sex', 'age', 'date_diff',
                 'head_size', 'head_location1', 'head_location2','head_location3','head_location4',
                 'qc_discrepancy', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5')

### Coefficient tests for each brain measurement and trait with Wald test and adj. p-values
#pvalues_lrtest_neuroticism <- list()
pvalues_ttest <- list()
# confounder matrix
Z <- as.matrix(mri_neuroticism_dummy[,c(confounders)])
# Loop over each brain structure
for (brain_structure in brain_structures) {
  print(brain_structure)
  #regress out head location in scanner
  for (trait in personality_traits) {
    X <- as.matrix(mri_neuroticism_dummy[brain_structure])
    Y <- as.matrix(mri_neuroticism_dummy[trait])
    # Wald test
    cit_params <- list(cit='wald')
    cit_params$'params_cit' <- list(lm_formula=stats::as.formula(paste(trait, "~ . + I(age^2)+age:sex+I(date_diff^2)")))
    wald_test <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                                   cit_with_parameters = cit_params)

    # Store the test result in the list
    pvalues_ttest[[paste0(brain_structure, trait)]] <- wald_test$p
    #pvalues_lrtest_neuroticism[[paste0(brain_structure, trait)]] <- pvalue_lrtest$`Pr(>Chisq)`[2]
  }
}
## p-values
p_val_df <- do.call(rbind, lapply(pvalues_ttest, as.data.frame))
p_val_df$p_adj <- stats::p.adjust(p_val_df[,1], method='BH')
colnames(p_val_df) <- c('p_unadj', 'p_adj')
p_val_df[order(p_val_df$p_unadj),]
#data.table::fwrite(p_val_df, file = path_to_pval_brain_trait_each)

### Manhattan plot
p_df_adj_plot <- p_val_df['p_unadj']
p_df_adj_plot <- p_df_adj_plot %>%
  tibble::rownames_to_column(var = "Label") %>%
  tidyr::extract(Label, into = c("Structure", "Trait"),
                 regex = "brain_structure_([0-9]+)([A-Za-z]+)")
#get ids of significant p-values from Avinun 2020
temporal_gyrus_dilligence <- paste0(brain_structures[which(ids_brain_structure=='26750')], 'diligence')
thalamus_dilligence <- paste0(brain_structures[which(ids_brain_structure=='26558')], 'diligence')
postcentral_gyrus_curiosity <- paste0(brain_structures[which(ids_brain_structure=='26776')], 'curiosity')
cerebral_peduncle_nervousness <- paste0(brain_structures[which(ids_brain_structure=='25071')], 'nervousness')
transverse_gyrus_dilligence <- paste0(brain_structures[which(ids_brain_structure=='26753')], 'diligence')
structure_trait_significant_avinun <- c(temporal_gyrus_dilligence, thalamus_dilligence, postcentral_gyrus_curiosity,
                                        cerebral_peduncle_nervousness, transverse_gyrus_dilligence)

# Convert Structure to numeric (important if not already numeric)
p_structure_trait <- p_df_adj_plot %>%
  mutate(Structure = as.numeric(Structure),
         Trait = as.factor(Trait)) %>%
  arrange(Trait, Structure) %>%
  mutate(X_Axis=1:nrow(p_df_adj_plot)) %>%

  # Add highlight and annotation information
  mutate( is_annotate=ifelse(-log10(p_unadj)>2.5, "yes", "no"))
p_structure_trait$combined <- with(p_structure_trait, paste("brain_structure_", Structure, Trait, sep = ""))
p_structure_trait$is_highlight <- ifelse(p_structure_trait$combined %in% structure_trait_significant_avinun, 'yes','no')

# labels for small p-values
ids_freesurfer_atlas <- data.table::fread(file=path_to_freesurfer_dk_atlas_ids,header=TRUE)
ids_freesurfer_aseg <- stats::na.omit(data.table::fread(file=path_to_freesurfer_aseg_ids,header=TRUE))
ids_FA_all <- stats::na.omit(data.table::fread(file=path_to_fractional_anisotropy_ids,header=TRUE))
ids_FA_lr <- ids_FA_all[grepl('left', ids_FA_all$Description) | grepl('right', ids_FA_all$Description) ,]
ids_surface_area_lr <- ids_freesurfer_atlas[1:68,]
ids_cortical_thick_lr <- ids_freesurfer_atlas[69:136,]
ids_gray_matter_vols_lr <- ids_freesurfer_aseg[c(47:50,61:64,67:68,72:73,79:82, 86:87),]
ids_lr <- rbind(ids_cortical_thick_lr,ids_surface_area_lr, ids_gray_matter_vols_lr, ids_FA_lr)
#ids_brain <- rbind(ids_lr[seq(from = 1, to = nrow(ids_lr), by = 2),], ids_whole_brain_aseg)
#lab_small_p_df <- ids_brain[p_structure_trait[which(p_structure_trait$is_annotate == 'yes'),]$Structure,]
labels <- c('medialorbitofrontal', 'rostralanteriorcingulate', 'posterior thalamic radiation',
            'isthmuscingulate', 'Caudate', 'medial lemniscus', 'posterior thalamic radiation',
            'rostralanteriorcingulate', 'pontine crossing tract')
#bonferroni adjusted significance level for 642 individual tests and significance level of 0.05 (on log-scale)
log_p_bonferroni_individual <- -log10(0.05/642)


# Plotting similar to a Manhattan plot
plot_each_struc_trait <- ggplot2::ggplot(p_structure_trait, ggplot2::aes(x = X_Axis, y = -log10(p_unadj), color = Trait), alpha=0.8, size=1.3) +
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::scale_x_continuous(labels = unique(p_structure_trait$Trait), breaks = seq(53, max(p_structure_trait$X_Axis), by = 107)) +
  ggplot2::scale_color_manual(values = rep(c(palet_discrete[7], palet_discrete[5]), 22 )) +
  ggplot2::labs(x = "Trait - Each brain structure",
       y = expression("-log"[10] * "(p)"),
       color = "Trait") +
  # Custom the theme:
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = element_text(size = 20),  # Change x-axis label size
    axis.title.y = element_text(size = 20),  # Change y-axis label size
    axis.text.x = element_text(size = 12),   # Change x-axis tick label size
    axis.text.y = element_text(size = 12),   # Change y-axis tick label size
    legend.position="none",
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank()
  ) +
  ggplot2::ylim(0,4.2) +
  # add colum separation for each trait
  geom_vline(xintercept = c(0,107, 214, 321, 428, 535, 642), linetype = "dashed", color = "black") +
  geom_hline(yintercept=c(0, log_p_bonferroni_individual), linetype=c("dashed", "solid"))+
  # Add highlighted points
  ggplot2::geom_point(data=subset(p_structure_trait, is_highlight=="yes"), color="black", size=2.5) +
  # Add highlighted points
  ggplot2::geom_point(data=subset(p_structure_trait, is_annotate=="yes"), color="red2", size=2.5) +
  # Add label using ggrepel to avoid overlapping
  ggrepel::geom_label_repel( data=subset(p_structure_trait, is_annotate=="yes"),  ggplot2::aes(label=labels), size=4)
#ggplot2::ggsave(paste0(path_to_save_plots, 'individual_p_values_wald.png'), plot_each_struc_trait, width = 7, height = 7, dpi = 300)

#### joint significance of brain structures
## Wald test for each trait separately
lm_models_joint <- list()
p_joint_wald <- list()
Z <- as.matrix(mri_neuroticism_dummy[,c(confounders)])
for (trait in personality_traits){
  X <- as.matrix(mri_neuroticism_dummy[brain_structures])
  Y <- as.matrix(mri_neuroticism_dummy[trait])
  # Wald test
  cit_params <- list(cit='wald')
  cit_params$'params_cit' <- list(lm_formula=stats::as.formula(paste(trait, "~ . + I(age^2)+age:sex+I(date_diff^2)")))
  wald_test <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                                 cit_with_parameters = cit_params)
  p_joint_wald[[trait]] <- wald_test$p
}
p_joint_df <- do.call(rbind, lapply(p_joint_wald, as.data.frame))
p_joint_df$p_adj <- stats::p.adjust(p_joint_df[,1], method='BH')
colnames(p_joint_df) <- c('p_unadj', 'p_adj')
p_total_wald <- rbind(p_val_df, p_joint_df)
p_total_wald$p_adj <- stats::p.adjust(p_total_wald[,1], method='BH')
#data.table::fwrite(p_total_wald, file = path_to_pval_brain_trait_each_w_joint)

### Deep-RCoT
rcots_p_joint_traits <- list()
Z <- sapply(as.data.frame(mri_neuroticism_dummy[, confounders]), as.numeric)
X <- as.matrix(mri_neuroticism_dummy[, brain_structures])
# for each trait
for (trait in personality_traits){
  print(trait)
  Y <- as.matrix(mri_neuroticism_dummy[, trait])
  rcots_p_joint_trait <- rep(1,20)
  for(seed in 1:20){
    cit_params <- list(cit='RCOT', params_cit=list(seed=seed))
    rcot_p <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                           cit_with_parameters = cit_params)$p
    rcots_p_joint_trait[seed] <- rcot_p
  }
  rcots_p_joint_traits[[trait]] <- mean(rcots_p_joint_trait)
}
# for all traits together
Y <- as.matrix(mri_neuroticism_dummy[, personality_traits])
rcots_p_joint_all_traits <- rep(1,20)
for(seed in 1:20){
  cit_params <- list(cit='RCOT', params_cit=list(seed=seed))
  rcot_p <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                         cit_with_parameters = cit_params)$p
  rcots_p_joint_all_traits[seed] <- rcot_p
}
rcots_p_joint_traits[['all_traits']] <- mean(rcots_p_joint_all_traits)
#data.table::fwrite(as.data.frame(rcots_p_joint_traits), file = path_to_pval_brain_trait_joint_rcot)

##plot joint p-values also as Manhattan plot
p_joint_df <- p_joint_df %>%
  tibble::rownames_to_column(var = "Trait") %>%
  dplyr::arrange(Trait) %>%
  dplyr::mutate(Test=rep('Wald', length(personality_traits))) %>%
  dplyr::mutate(p_adj=NULL)
p_joint_rcot <- t(as.data.frame(rcots_p_joint_traits))
p_joint_rcot <- data.frame(p_joint_rcot) %>%
  tibble::rownames_to_column(var = "Trait")%>%
  arrange(Trait) %>%
  dplyr::mutate(Test=rep('RCoT', length(rcots_p_joint_traits)))
colnames(p_joint_rcot) <- c('Trait', 'p_unadj', 'Test')
p_joint <- rbind(p_joint_df, p_joint_rcot)%>%
  # Add highlight and annotation information
  mutate( is_annotate=ifelse(-log10(p_unadj)>2.5, "yes", "no"))
p_joint$X_Axis <- 1:nrow(p_joint)

plot_joint <- ggplot2::ggplot(p_joint, ggplot2::aes(x = X_Axis, y = -log10(p_unadj), color = Test), alpha=0.8, size=1.3) +
  ggplot2::geom_point(alpha = 0.5) +
  ggplot2::scale_x_continuous(labels = unique(p_joint$Test), breaks = seq(3.5, max(p_joint$X_Axis), by = 7)) +
  ggplot2::scale_color_manual(values = rep(c(palet_discrete[8], palet_discrete[3]), 22 )) +
  ggplot2::labs(x = "Test-Trait-All brain structures",
                y = expression("-log"[10] * "(p)"),
                color = "Test") +
  # Custom the theme:
  ggplot2::theme_bw() +
  # add colum separation for each trait
  geom_vline(xintercept = c(0,6,13), linetype = "dashed", color = "black") +
  geom_hline(yintercept=0)+
  ggplot2::theme(
    axis.title.x = element_text(size = 20),  # Change x-axis label size
    axis.title.y = element_text(size = 20),  # Change y-axis label size
    axis.text.x = element_text(size = 15),   # Change x-axis tick label size
    axis.text.y = element_text(size = 15),   # Change y-axis tick label size
    legend.position="none",
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank()
  )+
  ggplot2::ylim(0,4)+
  ggplot2::geom_point(data=subset(p_joint, is_annotate=="yes"), color=palet_discrete[10], size=2) +
  # Add label using ggrepel to avoid overlapping
  ggrepel::geom_label_repel( data=subset(p_joint, is_annotate=="yes"),  ggplot2::aes(label=Trait), size=4)
#ggplot2::ggsave(paste0(path_to_save_plots, 'joint_p_values.png'), plot_joint, width = 7, height = 7, dpi = 300)


### Deep-RCoT with fastsurfer embedding map for all traits
#fastsurfer ids
ids_IDPs <- data.table::fread(path_to_fastsurfer_ids, select='Field ID')
ids_IDPs_full_str <- ids_IDPs %>%
  dplyr::mutate_all(list(~ stringr::str_c(., '-2.0'))) %>%
  dplyr::pull(`Field ID`)
ids_fastsurfer <- ids_IDPs_full_str[1:139]
brain_structures_fast <- paste0("brain_structure_", 1:length(ids_fastsurfer))
##Variables used and data subset
ids_mri_trait_deep_rcot <- c(ids_confounders, ids_fastsurfer, unique(ids_personality))
mrifast_trait <- stats::na.omit(subset_ukb_data_mri[,ids_mri_trait_deep_rcot])

###Outputs
# Neuroticism
mrifast_trait$neuroticism <- scale(mrifast_trait[, id_neuroticism])
# Sociability score
mrifast_trait$sociability <- scale(rowSums(mrifast_trait[, ids_social]))
# Warmth score
mrifast_trait$warmth <- scale(rowSums(mrifast_trait[, ids_warmth]))
# Diligence score
mrifast_trait$diligence <- scale(rowSums(mrifast_trait[, ids_diligence]))
# Curiosity score
mrifast_trait$curiosity <- scale(rowSums(mrifast_trait[, ids_curiosity]))
# Nervousness score
mrifast_trait$nervousness <- scale(rowSums(mrifast_trait[, ids_nervousness]))

## remove all original columns used to obtain the scores
mrifast_trait <- mrifast_trait[, !names(mrifast_trait) %in% unique(ids_personality)]
personality_traits <- c('neuroticism', 'sociability', 'warmth', 'diligence', 'curiosity', 'nervousness')
colnames(mrifast_trait) <- c('sex', 'age', 'site', 'date',
                               'head_size', 'head_location1', 'head_location2','head_location3','head_location4',
                               'qc_discrepancy', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5',
                             brain_structures_fast, personality_traits)
## dummy confounder site
dummy_df_tmp <- fastDummies::dummy_cols(mrifast_trait, select_columns = 'site')
mrifast_trait_dummy <- dummy_df_tmp[, !names(dummy_df_tmp) %in% c('site', 'site_10003')]
## date difference to first date
mrifast_trait_dummy$date_diff <- mrifast_trait_dummy[,'date'] - min(mrifast_trait_dummy[,'date'])
mrifast_trait_dummy <- mrifast_trait_dummy[, !names(mrifast_trait_dummy) %in% c('date')]
#data.table::fwrite(mrifast_trait_dummy, file = path_to_save_ukb_fast_behavior)

# confounder names
confounders <- c(grep("^site_", colnames(mrifast_trait_dummy), value=TRUE),
                 'sex', 'age', 'date_diff',
                 'head_size', 'head_location1', 'head_location2','head_location3','head_location4',
                 'qc_discrepancy', 'gene1', 'gene2', 'gene3', 'gene4', 'gene5')

### Deep-RCoT
rcots_p_joint_fast <- list()
Z <- sapply(as.data.frame(mrifast_trait_dummy[, confounders]), as.numeric)
X <- as.matrix(mrifast_trait_dummy[brain_structures_fast])
# for each trait
for (trait in personality_traits){
  print(trait)
  Y <- as.matrix(mrifast_trait_dummy[, trait])
  rcots_p_joint_fast_ <- rep(1,20)
  for(seed in 1:20){
    cit_params <- list(cit='RCOT', params_cit=list(seed=seed))
    rcot_p <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                           cit_with_parameters = cit_params)$p
    rcots_p_joint_fast_[seed] <- rcot_p
  }
  rcots_p_joint_fast[[trait]] <- mean(rcots_p_joint_fast_)
}
# for all traits together
Y <- as.matrix(mrifast_trait_dummy[, personality_traits])
rcots_p_joint_all_traits_fast <- rep(1,20)
for(seed in 1:20){
  cit_params <- list(cit='RCOT', params_cit=list(seed=seed))
  rcot_p <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                         cit_with_parameters = cit_params)$p
  rcots_p_joint_all_traits_fast[seed] <- rcot_p
}
rcots_p_joint_fast[['all_traits']] <- mean(rcots_p_joint_all_traits_fast)

#plot next to joint tests
p_joint_rcot_fast <- t(as.data.frame(rcots_p_joint_fast))
p_joint_rcot_fast <- data.frame(p_joint_rcot_fast) %>%
  tibble::rownames_to_column(var = "Trait")%>%
  arrange(Trait) %>%
  dplyr::mutate(Test=rep('Deep-RCoT', length(p_joint_rcot_fast)))
colnames(p_joint_rcot_fast) <- c('Trait', 'p_unadj', 'Test')
p_joint <- rbind(p_joint_df, p_joint_rcot, p_joint_rcot_fast)%>%
  # Add highlight and annotation information
  mutate( is_annotate=ifelse(-log10(p_unadj)>2.5, "yes", "no"))
p_joint$X_Axis <- 1:nrow(p_joint)
log_p_bonferroni_joint <- -log10(0.05/6)
p_joint <- p_joint %>% mutate(Test = recode(Test,
                                 "Wald"="Freesurfer-Wald",
                                 "RCoT"="Freesurfer-RCoT",
                                 "Deep-RCoT"="Fastsurfer-RCoT"))
p_joint <- p_joint %>% mutate(all_traits = ifelse(Trait == 'all_traits', 'yes', 'no'))

plot_joint <- ggplot2::ggplot(p_joint, ggplot2::aes(x = X_Axis, y = -log10(p_unadj), color = Test)) +
  ggplot2::geom_point(size=3, color=c(rep(palet_discrete[7], 6), rep(palet_discrete[1], 14)),
                      shape=c(rep(17, 6), 19, rep(17,6), 19, rep(15,6))) + #, rep(palet_discrete[3], 7))) +
  ggplot2::scale_x_continuous(labels = unique(p_joint$Test), breaks = seq(3.5, max(p_joint$X_Axis), by = 7)) +
  #ggplot2::scale_color_manual(values = c(rep(palet_discrete[9], 6), rep(palet_discrete[2], 7), rep(palet_discrete[3], 7))) +
  ggplot2::labs(x = "DNCITs - Trait(s) - All brain structures",
                y = expression("-log"[10] * "(p)"),
                color = "Test") +
  # Custom the theme:
  ggplot2::theme_bw() +
  # add colum separation for each trait
  geom_vline(xintercept = c(0.5,6.5,13.5,20.5), linetype = "dashed", color = "black") +
  geom_hline(yintercept=c(0, log_p_bonferroni_joint), linetype=c("dashed", "solid"))+
  ggplot2::theme(
    axis.title.x = element_text(size = 20),  # Change x-axis label size
    axis.title.y = element_text(size = 20),  # Change y-axis label size
    axis.text.x = element_text(size = 15),   # Change x-axis tick label size
    axis.text.y = element_text(size = 15),   # Change y-axis tick label size
    legend.position="none",
    panel.border = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank()
  )+
  ggplot2::ylim(0,4.2)+
  ggplot2::geom_point(data=subset(p_joint, all_traits=="yes"), color=c(palet_discrete[10], palet_discrete[10]), size=3.5, shape=19) +
  ggplot2::geom_point(data=subset(p_joint, is_annotate=="yes"), color="red2", size=3.5) +
  # Add label using ggrepel to avoid overlapping
  ggrepel::geom_label_repel( data=subset(p_joint, is_annotate=="yes"),  ggplot2::aes(label=Trait), size=4)
#ggplot2::ggsave(paste0(path_to_save_plots, 'deep_joint_p_values.png'), plot_joint, width = 7, height = 7, dpi = 300)



##### all p-values with adjustment (not used)
p_ind <- tibble::as_tibble(p_val_df['p_unadj'], rownames='Trait') %>%
  mutate(Test='Wald_individual')
p_all <- rbind(p_ind, p_joint[,1:3])
p_all$p_adj <- stats::p.adjust(p_all$p_unadj, method='BH')
p_all <- p_all %>%
  arrange(p_unadj)

## qq-plots to compare to uniform p-values
p_all$log_p_unadj = -log10(p_all$p_unadj)
tests <- c('Wald_individual', 'Wald', 'RCoT', 'Deep-RCoT')
qqplots <- list()
ks_tests <- list()
for(test in tests){
  p_all_test <- p_all[p_all$Test == test,]
  # Generate the expected theoretical quantiles for a uniform distribution
  n = nrow(p_all_test)
  theoretical_quantiles = (1:n - 0.5) / n
  theoretical_quantiles = sort(-log10(theoretical_quantiles))
  # Sort the observed -log10(p_unadj) values
  observed_quantiles = sort(p_all_test$log_p_unadj)
  # Create the QQ plot
  qqplots[[test]] <- ggplot2::ggplot(data = data.frame(x = theoretical_quantiles, y = observed_quantiles),
                           ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position="none",
      panel.border = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank()
    )+
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue") +
    ggplot2::xlab(expression("Theoretical Quantiles -log"[10] *"(Uniform))")) +
    ggplot2::ylab(expression("Observed Quantiles -log"[10] * "(p)"))
  ks_tests[[test]] <- ks.test(p_all_test$p_unadj, punif)
}
# Print the plot
#ggplot2::ggsave(paste0(path_to_save_plots, 'qqplot_ind_wald.png'), qqplots[['Wald_individual']], width = 7, height = 7, dpi = 300)

