##### We estimate the confounding effects of Z on
# 1) the brain structures obtained via the Fastsurfer embedding map
# 2) the behavioural traits
# using the Avinun 2020 UKB data
# and GAMs for each model with flexible smooth terms for continuous confounders
library(mgcv)

path_to_ukb_folder <- 'M:/CITs/Application/UKB_data/'
path_to_save_ukb_avinun <- paste0(path_to_ukb_folder, 'ukb_avinun.csv')
path_to_save_ids_brain_avinun <- paste0(path_to_ukb_folder, 'ids/ids_brain_avinun.csv')

# Avinun 2020 ukb data
ukb_avinun <- as.data.frame(data.table::fread(file=path_to_save_ukb_avinun,header=TRUE))

# List of confounders
confounders <- c("site_11001", "site_11002", "site_11003", "site_11004", "site_11005",
                 "site_11006", "site_11007", "site_11008", "site_11009", "site_11010",
                 "site_11011", "site_11012", "site_11013", "site_11014", "site_11016",
                 "site_11017", "site_11018", "site_11020", "site_11021", "site_11022",
                 "site_11023", "sex", "age", "date_diff", "head_size", "head_location1",
                 "head_location2", "head_location3", "head_location4", "qc_discrepancy",
                 "gene1", "gene2", "gene3", "gene4", "gene5")

# Binary confounders
binary_confounders <- c("sex", grep("^site", confounders, value = TRUE))

# Continuous confounders
continuous_confounders <- setdiff(confounders, binary_confounders)

# Prepare the list of brain structures
ids_brain_structure <- data.table::fread(file=path_to_save_ids_brain_avinun,header=TRUE)$ids_brain_structure
brain_structures <- paste0("brain_structure_", 1:length(ids_brain_structure))
#example with only 3 brain structures
brain_structures <- brain_structures[1:3]

# Fit a GAM for each brain structure and analyze the effects
results_brain_structures <- list()

for (brain_structure in brain_structures) {
  # Create the formula with smooth terms for continuous confounders
  smooth_terms <- paste("s(", continuous_confounders, ")", collapse = " + ")
  formula_str <- paste(brain_structure, "~", paste(binary_confounders, collapse = " + "), "+", smooth_terms)
  formula <- as.formula(formula_str)

  # Fit the GAM
  model <- gam(formula, data = ukb_avinun)

  # Store the summary of the model
  results_brain_structures[[brain_structure]] <- summary(model)
}

# Fit a GAM for each behavioural trait and analyze the effects
results_traits <- list()
personality_traits <- c('neuroticism', 'sociability', 'warmth', 'diligence', 'curiosity', 'nervousness')

for (trait in personality_traits) {
  # Create the formula with smooth terms for continuous confounders
  smooth_terms <- paste("s(", continuous_confounders, ")", collapse = " + ")
  formula_str <- paste(trait, "~", paste(binary_confounders, collapse = " + "), "+", smooth_terms)
  formula <- as.formula(formula_str)

  # Fit the GAM
  model <- gam(formula, data = ukb_avinun)

  # Store the summary of the model
  results_traits[[trait]] <- summary(model)
}

