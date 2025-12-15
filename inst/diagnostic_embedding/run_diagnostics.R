#### Diagnostic for different embeddings across CI conditions
#### Seeds 1:100, different eps_sigmaY configurations
# Load required packages and functions
library(dplyr)
library(devtools)
load_all()  # Load the dncitPaper package

# Configuration
idx_samples <- 1:4
n_sample = list(460, 1100, 5000, 10000)
xz_modes <- c('Sigma=I_p')
seeds <- 1:100
Y_age <- FALSE
standardize_ridge_lasso <- TRUE

# Define configurations: (ci_condition, eps_sigmaY, embedding_obs)
# eps_sigmaY = 0.5: both CI and No_CI with scratch and medicalnet_ft
# eps_sigmaY = 0.1, 1.0: No_CI only with scratch
configs <- list(
  list(ci_condition = '/CI/', eps_sigmaY = 0.5, embedding_obs = c('scratch', 'medicalnet_ft')),
  list(ci_condition = '/No_CI/', eps_sigmaY = 0.5, embedding_obs = c('scratch', 'medicalnet_ft')),
  list(ci_condition = '/No_CI/', eps_sigmaY = 0.1, embedding_obs = 'scratch'),
  list(ci_condition = '/No_CI/', eps_sigmaY = 1.0, embedding_obs = 'scratch')
)

# Run diagnostics for each configuration
for(config in configs){
  ci_condition <- config$ci_condition
  eps_sigmaY <- config$eps_sigmaY
  embedding_obs <- config$embedding_obs
  
  cat(sprintf("\n=== Running configuration: CI=%s, eps_sigmaY=%s, embeddings=%s ===\n", 
              ci_condition, eps_sigmaY, paste(embedding_obs, collapse=", ")))
  
  for(seed in seeds){
    cat(sprintf("Processing seed %d/%d\n", seed, max(seeds)))
    
    # Cache baseline results per seed (constant across sample sizes)
    baseline_cache <- NULL
    
    for(xz_mode in xz_modes){
      for(idx_sample in idx_samples){
        if(idx_sample == 1){
          # Compute baseline results for first sample size and cache them
          results <- auto_diagnostic(
            experiment_dir = paste0("/sc/home/marco.simnacher/ukbiobank/data", ci_condition, 
                                   n_sample[[idx_sample]], "/", seed, "/eps_sigmaY=", 
                                   as.character(eps_sigmaY)),
            embedding_obs = embedding_obs,
            seed = seed,
            debug_Y = FALSE,
            lambda_choice = "1se",
            Y_age = Y_age,
            standardize_ridge_lasso = standardize_ridge_lasso,
            xz_mode = xz_mode
          )
          baseline_cache <- results$baseline_results
        } else {
          # Reuse cached baseline results for subsequent sample sizes
          results <- auto_diagnostic(
            experiment_dir = paste0("/sc/home/marco.simnacher/ukbiobank/data", ci_condition, 
                                   n_sample[[idx_sample]], "/", seed, "/eps_sigmaY=", 
                                   as.character(eps_sigmaY)),
            embedding_obs = embedding_obs,
            seed = seed,
            debug_Y = FALSE,
            baseline_results_cached = baseline_cache,
            lambda_choice = "1se",
            Y_age = Y_age,
            standardize_ridge_lasso = standardize_ridge_lasso,
            xz_mode = xz_mode
          )
        }
      }
    }
  }
}

cat("\n=== All diagnostics completed ===\n")

