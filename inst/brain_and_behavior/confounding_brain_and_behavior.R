#### Confounder and DNCITs
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
path_to_fastsurfer_ids <- paste0(path_to_ukb,'IDPs/Ids_IDPs.csv')
path_to_fractional_anisotropy_ids <- paste0(path_to_ukb,'Real-world_application/ids_FA_measures.csv')
path_to_freesurfer_dk_atlas_ids <- paste0(path_to_ukb,'Real-world_application/ids_freesurfer_Desikan-Killiany_atlas.csv')
path_to_freesurfer_aseg_ids <- paste0(path_to_ukb,'Real-world_application/ids_freesurfer_ASEG.csv')
path_to_save_ids_brain_avinun <- paste0(path_to_ukb,'ids/ids_brain_avinun.csv')
path_to_save_ids_confounder_avinun <- paste0(path_to_ukb,'ids/ids_confounder_avinun.csv')
path_to_save_ids_personality_avinun <- paste0(path_to_ukb, 'ids/ids_personality_avinun.csv')
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
palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")

### Subset available UKB to obtain available col names
ukb_whole_columns <- data.table::fread(file=path_to_ukb_data, header=TRUE, nrows=10)

## load preprocessed data and ids
ids_brain_structure <- data.table::fread(file=path_to_save_ids_brain_avinun,header=TRUE)$ids_brain_structure
brain_structures <- paste0("brain_structure_", 1:length(ids_brain_structure))
ids_confounders <- data.table::fread(file=path_to_save_ids_confounder_avinun,header=TRUE)$ids_confounder
ids_personality <- data.table::fread(file=path_to_save_ids_personality,header=TRUE)$ids_personality
ukb_avinun <- as.data.frame(data.table::fread(file=path_to_save_ukb_avinun,header=TRUE))
ids_all <- c(ids_confounders, ids_brain_structure, ids_personality)


### 1.1) Confounder control: test if nonlinear relationships with confounders are satisfactory addressed in Wald tests (Avinun+RCoT)
confounders_short <- c('age', 'sex', 'genes', 'head_size', 'head_position', 'date_diff', 'site')
brain_structures_lin_unconfounded <- list()
# for each confounder compute linearly unconfounded brain structure measurements
for(confounder in confounders_short){
  count <- 1
  brain_structures_lin_unconfounded[[confounder]] <- matrix(nrow = nrow(ukb_avinun), ncol=length(brain_structures))
  for(brain_structure in brain_structures){
    if (confounder == 'head_position'){
      lm_formula <-stats::as.formula(paste(brain_structure,'~head_location1+head_location2+head_location3+head_location4'))
    }else if (confounder == 'genes'){
      lm_formula <-stats::as.formula(paste(brain_structure,'~gene1+gene2+gene3+gene4+gene5'))
    }else if (confounder == 'site'){
      lm_formula <-stats::as.formula(paste(brain_structure,'~', paste(grep("^site_", colnames(ukb_avinun), value=TRUE), collapse = "+")))
    }else if(confounder == 'age'){
      lm_formula <-stats::as.formula(paste(brain_structure,'~', confounder, '+I(age^2)+age:sex'))
    }else if(confounder == 'date_diff'){
      lm_formula <-stats::as.formula(paste(brain_structure,'~', confounder, '+I(date_diff^2)'))
    }else if(confounder == 'sex'){
      lm_formula <-stats::as.formula(paste(brain_structure,'~', confounder, '+age:sex'))
    }else{
      lm_formula <-stats::as.formula(paste(brain_structure,'~', confounder))
    }

    lm_model <- lm(lm_formula, ukb_avinun)
    brain_structures_lin_unconfounded[[confounder]][, count] <- lm_model$residuals
    count <- count +1
  }
  colnames(brain_structures_lin_unconfounded[[confounder]]) <- brain_structures
}
# CIT unconfounded residuals \indep confounder | other_confounders
mri_lin_unconf <- list()
rcots_p_mean <- list()
#wald_p_confounder <- list()
for (confounder in confounders_short){
  print(confounder)
  mri_lin_unconf[[confounder]] <- cbind(ukb_avinun[, !names(ukb_avinun) %in% brain_structures], brain_structures_lin_unconfounded[[confounder]])
  X <- as.matrix(brain_structures_lin_unconfounded[[confounder]])
  if (confounder == 'head_position'){
    head_pos <- c('head_location1','head_location2','head_location3','head_location4')
    Z <- sapply(as.data.frame(ukb_avinun[confounders[confounders %notin% head_pos]]), as.numeric)
    Y <- as.matrix(ukb_avinun[, head_pos])
  }else if (confounder == 'genes'){
    genes_ <- c('gene1','gene2','gene3','gene4','gene5')
    Z <- sapply(as.data.frame(ukb_avinun[confounders[confounders %notin% genes_]]), as.numeric)
    Y <- as.matrix(ukb_avinun[, genes_])
  }else if (confounder == 'site'){
    sites_ <- grep("^site_", colnames(ukb_avinun), value=TRUE)
    Z <- sapply(as.data.frame(ukb_avinun[confounders[confounders %notin% sites_]]), as.numeric)
    Y <- as.matrix(ukb_avinun[, sites_])
  }else{
    Z <- sapply(as.data.frame(ukb_avinun[confounders[confounders %notin% confounder]]), as.numeric)
    Y <- as.matrix(ukb_avinun[confounder])
  }
  rcots_p <- rep(1, 20)
  for(seed in 1:20){
    cit_params <- list(cit='RCOT', params_cit=list(seed=seed))
    rcots_p[seed] <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                                  cit_with_parameters = cit_params)$p
  }
  rcots_p_mean[[confounder]] <- mean(rcots_p)

  # cit_params <- list(cit='wald')
  # cit_params$'params_cit' <- list(lm_formula=stats::as.formula(paste(confounder, "~ . + I(age^2)+age:sex+I(date_diff^2)")))
  # wald_p_confounder[[confounder]] <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
  #                           cit_with_parameters = cit_params)$p
}
# still evidence for nonlinear confounding for age, genes, head size, head position (and potentially date diff)
#data.table::fwrite(as.data.frame(rcots_p_mean), file = path_to_pval_brain_trait_confounder_control)
to_latex_text <- function(lst) {
  latex_str <- ""
  for (name in names(lst)) {
    value <- lst[[name]]
    if (is.numeric(value)) {
      rounded_value <- signif(value, digits = 3)  # Rounding to 5 significant digits
      latex_str <- paste(latex_str, sprintf("%s %s", name, rounded_value), sep="; ")
    } else {
      latex_str <- paste(latex_str, sprintf("%s %s", name, value), sep="; ")
    }
  }
  return(latex_str)
}
to_latex_text(rcots_p_mean)


### 2) Additional confounders
X <- as.matrix(mri_neuroticism_dummy[, brain_structures])
confounder <- as.matrix(mri_neuroticism_dummy[, personality_traits])
Z <- sapply(as.data.frame(mri_neuroticism_dummy[, names(mri_neuroticism_dummy) %in% confounders[confounders!='age']]), as.numeric)
age <- mri_neuroticism_dummy[,'age']
