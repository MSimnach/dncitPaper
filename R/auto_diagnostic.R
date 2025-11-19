#' Automatic diagnostic for embeddings
#'
#' Evaluates trained and baseline embeddings on a separate diagnostic dataset,
#' regenerating Y with the same parameters as the original experiment.
#'
#' @param experiment_dir Path to experiment directory (e.g., "/No_CI/550/1/")
#' @param embedding_obs Vector of trained embedding names (e.g., c("scratch"), c("medicalnet_ft", "medicalnet_ft_frozen"))
#' @param diagnostic_csv Path to diagnostic dataset CSV. If NULL, defaults to "{UKB_PATH}/t1_paths_diag.csv"
#' @param seed Random seed for Y generation and glmnet (default: 1)
#' @param extract_trained Extract embeddings from checkpoint (default: TRUE)
#' @param baseline_embeddings Vector of baseline embeddings to evaluate (default: c("fastsurfer", "condVAE", "freesurfer", "medicalnet"))
#' @param alpha Elastic net mixing parameter (default: 1 for lasso)
#' @param test_prop Test set proportion (default: 0.2)
#' @param nfolds Number of CV folds (default: 10)
#' @param lambda_choice Lambda selection method: "min" or "1se" (default: "1se")
#' @param baseline_results_cached Pre-computed baseline results to reuse (default: NULL)
#' @param Y_age Generate Y based on age (default: FALSE)
#' @param xz_mode XZ mode: "Sigma=I_p", "independent", "correlated" (default: "Sigma=I_p")
#'
#' @return List containing results data frame, baseline results, and paths to saved files
#' @export
#'
#' @examples
#' \dontrun{
#' # Evaluate single trained embedding
#' results <- auto_diagnostic(
#'   experiment_dir = "/path/to/No_CI/550/1",
#'   embedding_obs = c("scratch")
#' )
#'
#' # Evaluate multiple trained embeddings
#' results <- auto_diagnostic(
#'   experiment_dir = "/path/to/No_CI/550/1",
#'   embedding_obs = c("scratch", "medicalnet_ft", "medicalnet_ft_frozen")
#' )
#' }
auto_diagnostic <- function(
  experiment_dir,
  embedding_obs,
  diagnostic_csv = NULL,
  seed = 1,
  extract_trained = TRUE,
  baseline_embeddings = c("fastsurfer", "condVAE", "freesurfer", "medicalnet"),#"pooled_brainsynth", "tucker_brainsynth"),
  alpha = 0.3,
  test_prop = 0.2,
  nfolds = 10,
  lambda_choice = "1se",
  standardize_ridge_lasso = TRUE,
  debug_Y = TRUE,
  baseline_results_cached = NULL,
  Y_age = FALSE,
  xz_mode = "Sigma=I_p"
) {
  
  # Set single thread for reproducibility
  data.table::setDTthreads(1)
  
  # ============================================================================
  # 1. Load Configuration and Setup
  # ============================================================================
  cat("\n=== Auto Diagnostic for Embeddings ===\n")
  cat("Experiment directory:", experiment_dir, "\n")
  cat("Trained embeddings:", paste(embedding_obs, collapse=", "), "\n\n")
  
  # Get UKB_PATH from environment
  path_to_ukb_data <- Sys.getenv("UKB_PATH", unset = NA)
  if (is.na(path_to_ukb_data) || path_to_ukb_data == "") {
    stop("Environment variable UKB_PATH is not set. Please define it in your .Renviron.")
  }
  
  # Set default diagnostic_csv if not provided
  if (is.null(diagnostic_csv)) {
    diagnostic_csv <- paste0(path_to_ukb_data, "/t1_paths_diag.csv")
  }
  
  existing_embeddings <- c()
  for (emb in embedding_obs) {
    emb_dir <- file.path(experiment_dir, emb)
    if (dir.exists(emb_dir)) {
      existing_embeddings <- c(existing_embeddings, emb)
    } else {
      warning("Trained embedding directory does not exist, will skip: ", emb_dir)
    }
  }
  
  # If no trained embeddings exist, skip the whole auto_diagnostic
  if (length(existing_embeddings) == 0) {
    cat("No trained embedding directories found. Skipping auto_diagnostic.\n")
    return(invisible(NULL))
  }
  
  cat("Found", length(existing_embeddings), "out of", length(embedding_obs), 
      "trained embedding directories\n")
  
  # Read config from first existing trained embedding
  first_existing_emb <- existing_embeddings[1]
  first_emb_dir <- file.path(experiment_dir, first_existing_emb)
  config_path <- file.path(first_emb_dir, "config.json")
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  
  config <- jsonlite::read_json(config_path)
  cat("Loaded configuration from:", config_path, "\n")
  
  # Extract parameters
  post_non_lin <- config$post_non_lin
  eps_sigmaX <- config$eps_sigmaX
  eps_sigmaY <- config$eps_sigmaY
  eps_sigmaZ <- config$eps_sigmaZ
  embedding_orig <- config$embedding_orig
  confounder <- config$confounder
  g_z <- config$g_z
  if (Y_age) {
    is_ci <- "Y_age"
  } else {
    is_ci <- config$is_ci
  }
  
  cat("Configuration: post_non_lin=", post_non_lin, ", eps_sigmaY=", eps_sigmaY, 
      ", confounder=", confounder, ", g_z=", g_z, ", xz_mode=", xz_mode, "\n", sep="")
  cat("\n")
  
  # Create output directory
  out_dir_diagnostic <- file.path(experiment_dir, "diagnostic")
  dir.create(out_dir_diagnostic, recursive = TRUE, showWarnings = FALSE)
  cat("Output directory:", out_dir_diagnostic, "\n\n")
  
  # Set seed
  set.seed(seed)
  
  # ============================================================================
  # 2. Load Diagnostic Dataset and Baseline Embeddings
  # ============================================================================
  cat("=== Loading Diagnostic Dataset and Baseline Embeddings ===\n")
  
  # Load diagnostic T1 paths
  if (!file.exists(diagnostic_csv)) {
    stop("Diagnostic CSV not found: ", diagnostic_csv)
  }
  t1_paths_diag <- data.table::fread(diagnostic_csv, header=TRUE, nThread = 1)
  cat("Loaded", nrow(t1_paths_diag), "samples from diagnostic dataset\n")
  
  # Initialize list to store all embeddings for ID matching
  all_datasets <- list(t1_paths_diag = t1_paths_diag)
  
  # Load baseline embeddings
  baseline_emb_list <- list()
  
  for (emb_name in baseline_embeddings) {
    cat("Loading baseline embedding:", emb_name, "... ")
    tryCatch({
      if (emb_name == "fastsurfer") {
        emb_data <- data.table::fread(paste0(path_to_ukb_data, "/ukb_fastsurfer.csv"), 
                                     header=TRUE, nThread = 1)
        colnames(emb_data)[1] <- "id"
        baseline_emb_list[[emb_name]] <- emb_data
        all_datasets[[emb_name]] <- emb_data
        cat("OK (", ncol(emb_data)-1, " features)\n", sep="")
        
      } else if (emb_name == "condVAE") {
        emb_data <- data.table::fread(paste0(path_to_ukb_data, "/ukb_condVAE.csv"), 
                                     header=FALSE, nThread = 1)
        colnames(emb_data)[1] <- "id"
        baseline_emb_list[[emb_name]] <- emb_data
        all_datasets[[emb_name]] <- emb_data
        cat("OK (", ncol(emb_data)-1, " features)\n", sep="")
        
      } else if (emb_name == "freesurfer") {
        emb_data <- data.table::fread(paste0(path_to_ukb_data, "/ukb_freesurfer.csv"), 
                                     header=TRUE, nThread = 1)
        colnames(emb_data)[1] <- "id"
        baseline_emb_list[[emb_name]] <- emb_data
        all_datasets[[emb_name]] <- emb_data
        cat("OK (", ncol(emb_data)-1, " features)\n", sep="")
        
      } else if (emb_name == "medicalnet") {
        medicalnet_dir <- paste0(path_to_ukb_data, "/medicalnet_embeddings")
        medicalnet_emb <- as.matrix(arrow::read_parquet(paste0(medicalnet_dir, "/embeddings.parquet")))
        medicalnet_idx <- data.table::fread(paste0(medicalnet_dir, "/embeddings_index.csv"), nThread = 1)
        emb_data <- cbind(id = medicalnet_idx$subject_id, as.data.frame(medicalnet_emb))
        baseline_emb_list[[emb_name]] <- emb_data
        all_datasets[[emb_name]] <- emb_data
        cat("OK (", ncol(emb_data)-1, " features)\n", sep="")
      } else if (emb_name == "boc_brainsynth") {
        emb_data <- arrow::read_parquet(paste0(path_to_ukb_data, '/brainsynth_embeddings/baseline_vqvae/embeddings/all_boc_embedding.parquet'))
        emb_data <- as.data.frame(emb_data)
        boc_idx <- data.table::fread(paste0(path_to_ukb_data, "/brainsynth_embeddings/baseline_vqvae/embeddings/subject_ids.csv"))
        emb_data <- cbind(id = boc_idx$x, as.data.frame(emb_data))
        baseline_emb_list[[emb_name]] <- emb_data
        all_datasets[[emb_name]] <- emb_data
        cat("OK (", ncol(emb_data)-1, " features)\n", sep="")
      } else if (emb_name == "pooled_brainsynth") {
        emb_data <- arrow::read_parquet(paste0(path_to_ukb_data, '/brainsynth_embeddings/baseline_vqvae/embeddings/all_pooled_embeddings.parquet'))
        pooled_idx <- data.table::fread(paste0(path_to_ukb_data, "/brainsynth_embeddings/baseline_vqvae/embeddings/subject_ids.csv"))
        emb_data <- as.data.frame(emb_data)
        emb_data <- cbind(id = pooled_idx$x, as.data.frame(emb_data))
        baseline_emb_list[[emb_name]] <- emb_data
        all_datasets[[emb_name]] <- emb_data
        cat("OK (", ncol(emb_data)-1, " features)\n", sep="")
      } else if (emb_name == "tucker_brainsynth") {
        emb_data <- arrow::read_parquet(paste0(path_to_ukb_data, '/brainsynth_embeddings/baseline_vqvae/embeddings/all_tucker_embeddings.parquet'))
        colnames(emb_data)[1] <- "id"
        emb_data <- as.data.frame(emb_data)
        baseline_emb_list[[emb_name]] <- emb_data
        all_datasets[[emb_name]] <- emb_data
        cat("OK (", ncol(emb_data)-1, " features)\n", sep="")
      }
    }, error = function(e) {
      cat("SKIP (not found)\n")
    })
  }
  
  # Find common IDs across all datasets
  common_ids <- Reduce(intersect, lapply(all_datasets, function(x) x$id))
  common_ids <- as.integer(common_ids) 
  if (length(common_ids) == 0) {
    stop("No common IDs found across datasets")
  }
  
  common_ids <- sort(common_ids)
  
  # Filter and reorder each dataset
  t1_paths_diag <- t1_paths_diag[match(common_ids, t1_paths_diag$id), ]
  
  for (emb_name in names(baseline_emb_list)) {
    baseline_emb_list[[emb_name]] <- baseline_emb_list[[emb_name]][
      match(common_ids, baseline_emb_list[[emb_name]]$id), 
    ]
  }

  # ADD VERIFICATION AFTER FILTERING (after line 233):
  cat("\n=== Verifying ID Matching After Filtering/Reordering ===\n")
  for (emb_name in names(baseline_emb_list)) {
    emb_ids <- baseline_emb_list[[emb_name]]$id
    
    # Check for NAs (shouldn't happen, but verify)
    if (any(is.na(emb_ids))) {
      warning("Found NA IDs in ", emb_name, " after filtering!")
    }
    
    # Check if IDs match common_ids (should be TRUE after reordering)
    ids_match <- all(emb_ids == common_ids, na.rm=TRUE)
    cat("  ", emb_name, ": ", nrow(baseline_emb_list[[emb_name]]), " rows, IDs match: ", ids_match, "\n", sep="")
    
    if (!ids_match) {
      # Show first few mismatches for debugging
      mismatch_idx <- which(emb_ids != common_ids)[1:min(5, length(which(emb_ids != common_ids)))]
      cat("    First mismatches at positions: ", paste(mismatch_idx, collapse=", "), "\n", sep="")
      cat("    Expected: ", paste(common_ids[mismatch_idx], collapse=", "), "\n", sep="")
      cat("    Got:      ", paste(emb_ids[mismatch_idx], collapse=", "), "\n", sep="")
    }
  }
  cat("\n")
  
  cat("All datasets filtered and reordered to match common IDs\n\n")
  
  # ============================================================================
  # 3. Extract Trained Embeddings (if extract_trained=TRUE)
  # ============================================================================
  trained_emb_list <- list()
  
  if (extract_trained) {
    cat("=== Extracting Trained Embeddings ===\n")
    
    # Save filtered diagnostic paths
    t1_paths_diag_file <- file.path(out_dir_diagnostic, "t1_paths_diag.csv")
    write.csv(t1_paths_diag, t1_paths_diag_file, row.names = FALSE)
    cat("Saved filtered diagnostic paths:", t1_paths_diag_file, "\n\n")
    
    # Loop through each trained embedding
    for (emb in embedding_obs) {
      # Check if embedding directory exists
      emb_dir <- file.path(experiment_dir, emb)
      if (!dir.exists(emb_dir)) {
        warning("Skipping extraction for ", emb, ": directory does not exist: ", emb_dir)
        next
      }

      if (file.exists(file.path(out_dir_diagnostic, emb, "embeddings.parquet"))) {
        cat("Embeddings already extracted for:", emb, "\n")
      } else {
        cat("Extracting embeddings for:", emb, "\n")
        
        checkpoint_path <- file.path(experiment_dir, emb, "best_model.ckpt")
        if (!file.exists(checkpoint_path)) {
          warning("Checkpoint not found for ", emb, ": ", checkpoint_path)
          next
        }
        
        # Call Python extraction script
        system2("python", c(
          "inst/learn_embedding/extract_embeddings.py",
          "--checkpoint", checkpoint_path,
          "--input_csv", t1_paths_diag_file,
          "--output_dir", out_dir_diagnostic,
          "--config_name", emb,
          "--batch_size", "16",
          "--num_workers", "4",
          "--amp"
        ))
      }
      
      # Load extracted embeddings
      emb_path <- file.path(out_dir_diagnostic, emb, "embeddings.parquet")
      idx_path <- file.path(out_dir_diagnostic, emb, "embeddings_index.csv")
      
      if (file.exists(emb_path) && file.exists(idx_path)) {
        trained_emb <- as.matrix(arrow::read_parquet(emb_path))
        trained_idx <- data.table::fread(idx_path, nThread = 1)
        trained_embedding <- cbind(id = trained_idx$subject_id, as.data.frame(trained_emb))
        
        # Reorder to match common IDs
        trained_embedding <- trained_embedding[match(common_ids, trained_embedding$id), ]
        trained_emb_list[[emb]] <- trained_embedding
        
        cat("  Loaded embeddings: ", nrow(trained_embedding), " samples, ", 
            ncol(trained_embedding)-1, " features\n\n", sep="")
      } else {
        warning("Failed to extract embeddings for ", emb)
      }
    }
  } else {
    cat("=== Skipping Trained Embeddings Extraction ===\n")
    cat("Loading pre-extracted embeddings...\n\n")
    
    # Load pre-existing embeddings
    for (emb in embedding_obs) {
      # Check if embedding directory exists
      emb_dir <- file.path(experiment_dir, emb)
      if (!dir.exists(emb_dir)) {
        warning("Skipping loading for ", emb, ": directory does not exist: ", emb_dir)
        next
      }
      
      emb_path <- file.path(out_dir_diagnostic, emb, "embeddings.parquet")
      idx_path <- file.path(out_dir_diagnostic, emb, "embeddings_index.csv")
      
      if (file.exists(emb_path) && file.exists(idx_path)) {
        trained_emb <- as.matrix(arrow::read_parquet(emb_path))
        trained_idx <- data.table::fread(idx_path, nThread = 1)
        trained_embedding <- cbind(id = trained_idx$subject_id, as.data.frame(trained_emb))
        
        # Reorder to match common IDs
        trained_embedding <- trained_embedding[match(common_ids, trained_embedding$id), ]
        trained_emb_list[[emb]] <- trained_embedding
        cat("  Loaded", emb, ":", nrow(trained_embedding), "samples,", 
            ncol(trained_embedding)-1, "features\n")
      } else {
        warning("Pre-extracted embeddings not found for ", emb)
      }
    }
    cat("\n")
  }
  
  # ============================================================================
  # 4. Load and Process Confounders
  # ============================================================================
  cat("=== Loading and Processing Confounders ===\n")
  
  Z <- load_Z(path_to_ukb_data, confounder)
  cat("Loaded confounders:", confounder, "(", ncol(Z)-1, "variables )\n")
  
  # Add noise
  epsZ <- matrix(stats::rnorm((nrow(Z)*(ncol(Z)-1)), 0, eps_sigmaZ), nrow=nrow(Z), ncol=ncol(Z)-1)
  Z[,c(-1)] <- Z[,c(-1)]+epsZ
  
  # Filter to common IDs
  Z <- Z[match(common_ids, Z$id), ]
  
  # Remove zero columns and multicollinearity in one-hot-encoding of sites
  site_columns <- grep("^site", colnames(Z), value = TRUE)
  for (site in site_columns) {
    site_sum <- sum(Z[[site]])
    if (site_sum == 0) {
      Z <- Z[, !names(Z) %in% site]
    }
  }
  
  site_columns <- grep("^site", colnames(Z), value = TRUE)
  site_sum <- sum(as.data.frame(Z)[site_columns])
  
  # Remove one site column to avoid multicollinearity
  if (site_sum == nrow(Z)) {
    for (site in site_columns) {
      if (is_binary(Z[[site]])) {
        Z <- Z[, !names(Z) %in% site]
        break
      }
    }
  }
  
  # Standardize continuous confounders
  Z[,c(-1)] <- Z[,c(-1)] %>% 
    dplyr::mutate(dplyr::across(dplyr::where(function(x) !is_binary(x)), scale))
  
  cat("Processed confounders:", ncol(Z)-1, "variables after preprocessing\n\n")
  
  # ============================================================================
  # 5. Generate Y on Diagnostic Data
  # ============================================================================
  cat("=== Generating Y on Diagnostic Data ===\n")
  
  # Get fastsurfer embeddings for Y generation (X_orig)
  if ("fastsurfer" %in% names(baseline_emb_list)) {
    fastsurfer_emb <- baseline_emb_list[["fastsurfer"]]
  } else {
    fastsurfer_emb <- data.table::fread(paste0(path_to_ukb_data, "/ukb_fastsurfer.csv"), 
                                       header=TRUE, nThread = 1)
    colnames(fastsurfer_emb)[1] <- "id"
    fastsurfer_emb <- fastsurfer_emb[match(common_ids, fastsurfer_emb$id), ]
  }
  
  # Scale fastsurfer embeddings
  fastsurfer_emb <- as.data.frame(fastsurfer_emb)
  fastsurfer_emb[,c(-1)] <- scale(fastsurfer_emb[,c(-1)])
  
  # Generate Y based on CI or No_CI
  if (is_ci == "CI" || is.null(is_ci)) {
    cat("Generating Y under Conditional Independence (CI)\n")
    Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z, xz_mode=xz_mode)
    Y_info <- NULL  # Add this line
  } else if (is_ci == "No_CI") {
    cat("Generating Y under No Conditional Independence (No_CI)\n")
    if (debug_Y) {
      Y_list <- y_from_xz(Z[,c(-1)], eps_sigmaY, X=as.matrix(fastsurfer_emb[,c(-1)]), 
                    gamma=0.5, post_non_lin=post_non_lin, g_z=g_z, debug=debug_Y, xz_mode=xz_mode)
      Y <- Y_list$Y
      Y_info <- Y_list$info
    } else {
      Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, X=as.matrix(fastsurfer_emb[,c(-1)]), 
                    gamma=0.5, post_non_lin=post_non_lin, g_z=g_z, xz_mode=xz_mode)
      Y_info <- NULL  # Add this line
    }
  } else if (is_ci == "Y_age") {
    if ("age" %in% colnames(Z)) {
      Y <- Z$age
      Y_info <- NULL
    } else {
      stop("Age column not found in Z")
    }
  }
  
  # Save Y
  Y_id <- data.frame(id = fastsurfer_emb$id, Y = Y)
  if (is_ci == "Y_age") {
    Y_path <- file.path(out_dir_diagnostic, 'Y_diag_age.csv')
  } else {
    Y_path <- file.path(out_dir_diagnostic, 'Y_diag.csv')
  }
  write.csv(Y_id, Y_path, row.names = FALSE)
  cat("Saved Y to:", Y_path, "\n")
  cat("Y statistics: mean=", round(mean(Y), 4), ", sd=", round(sd(Y), 4), "\n\n", sep="")
  # Save Y_info if available
  if (debug_Y) {
    # Option 1: Save as RDS (preserves R object structure)
    Y_info_path <- file.path(out_dir_diagnostic, 'Y_info.rds')
    saveRDS(Y_info, Y_info_path)
    cat("Saved Y_info to:", Y_info_path, "\n")
    
    # Option 2: Also save as JSON for readability
    Y_info_json_path <- file.path(out_dir_diagnostic, 'Y_info.json')
    jsonlite::write_json(Y_info, Y_info_json_path, pretty = TRUE, auto_unbox = TRUE)
    cat("Saved Y_info (JSON) to:", Y_info_json_path, "\n")
    
    # Print some info
    cat("Y_info contents:\n")
    cat("  - Correlation pX vs pZ:", round(Y_info$corr_pX_pZ, 4), "\n")
    cat("  - Mode:", Y_info$mode, "\n")
    cat("  - Gamma:", Y_info$gamma, "\n")
    cat("  - Tau:", Y_info$tau, "\n\n")
  }
  
  # ============================================================================
  # 6. Evaluate All Embeddings
  # ============================================================================
  cat("=== Evaluating Embeddings with Glmnet ===\n")
  
  results_list <- list()
  
  # Evaluate baseline embeddings (or use cached results)
  if (!is.null(baseline_results_cached)) {
    cat("\n--- Using Cached Baseline Embeddings ---\n")
    # Use pre-computed baseline results
    for (i in seq_len(nrow(baseline_results_cached))) {
      results_list[[length(results_list) + 1]] <- baseline_results_cached[i, ]
      cat("Using cached", baseline_results_cached$embedding_name[i], 
          ": R²=", round(baseline_results_cached$r2_test[i], 4), 
          ", MSE=", round(baseline_results_cached$mse_test[i], 4), "\n", sep="")
    }
  } else if (length(baseline_emb_list) > 0) {
    cat("\n--- Baseline Embeddings ---\n")
    for (emb_name in names(baseline_emb_list)) {
      cat("Evaluating", emb_name, "... ")
      
      emb_data <- baseline_emb_list[[emb_name]]
      # Scale embeddings (remove id column)
      X_raw <- scale(as.matrix(emb_data[,c(-1)]))
      # DEBUG: Check for zero variance BEFORE scaling
      X_before_scale <- as.matrix(emb_data[,c(-1)])
      zero_var_cols_before <- apply(X_before_scale, 2, function(x) {
        var_x <- var(x, na.rm=TRUE)
        is.na(var_x) || var_x == 0 || var_x < 1e-10
      })
      if (any(zero_var_cols_before)) {
        cat(sprintf("[DEBUG] Found %d zero/low-variance columns BEFORE scaling in %s: %s\n", 
                    sum(zero_var_cols_before),
                    emb_name,
                    paste(which(zero_var_cols_before), collapse=", ")))
      }
      
      # Remove columns with NA values (zero variance columns)
      na_cols <- colSums(is.na(X_raw)) > 0
      if (any(na_cols)) {
        na_col_names <- names(X_raw)[na_cols]
        cat(sprintf("[INFO] Removing %d NA columns from X_obs: %s\n", 
                  sum(na_cols), 
                  paste(na_col_names, collapse=", ")))
        X <- X_raw[, c(!na_cols)]  # Keep id column (TRUE) and non-NA columns
      }else{
        X <- X_raw
      }

      # DEBUG: Check for zero variance AFTER removal
      zero_var_cols_after <- apply(X, 2, function(x) {
        var_x <- var(x, na.rm=TRUE)
        is.na(var_x) || var_x == 0 || var_x < 1e-10
      })
      if (any(zero_var_cols_after)) {
        cat(sprintf("[WARNING] Still found %d zero/low-variance columns AFTER removal in %s: %s\n", 
                    sum(zero_var_cols_after),
                    emb_name,
                    paste(which(zero_var_cols_after), collapse=", ")))
      }
      
      # Add diagnostic output
      cat("\n  X dimensions: ", nrow(X), "x", ncol(X), 
          ", mean:", round(mean(X), 4), 
          ", sd:", round(sd(X), 4),
          ", range: [", round(min(X), 4), ",", round(max(X), 4), "]\n", sep="")
      cat("  First few values: ", paste(round(X[1, 1:min(5, ncol(X))], 4), collapse=" "), "\n")
      # DEBUG: Check Y variance
      Y_var <- var(Y, na.rm=TRUE)
      cat("  Y variance:", round(Y_var, 6), "\n")
      if (Y_var == 0 || is.na(Y_var)) {
        cat("  [WARNING] Y has zero variance!\n")
      }

      # Check correlation between features and Y (with error handling)
      cor_with_y <- apply(X, 2, function(x) {
        tryCatch({
          cor(x, Y, use="complete.obs")
        }, warning = function(w) {
          cat(sprintf("    [WARNING in cor()] Column %d: %s\n", 
                      which(apply(X, 2, identical, x))[1], 
                      conditionMessage(w)))
          return(NA)
        }, error = function(e) {
          cat(sprintf("    [ERROR in cor()] Column %d: %s\n", 
                      which(apply(X, 2, identical, x))[1], 
                      conditionMessage(e)))
          return(NA)
        })
      })
      cat("  Max |correlation| with Y:", max(abs(cor_with_y), na.rm=TRUE), "\n")
      cat("  Features with |cor| > 0.1:", sum(abs(cor_with_y) > 0.1, na.rm=TRUE), "\n")
      cat("  Features with NA correlation:", sum(is.na(cor_with_y)), "\n")
      # Run glmnet
      result <- ridge_lasso_glmnet(
        X = X, 
        y = Y,
        alpha = alpha,
        test_prop = test_prop,
        nfolds = nfolds,
        lambda_choice = lambda_choice,
        seed = seed, 
        standardize = FALSE # already happens above
      )
      cat("  Lambda used:", result$lambda, "\n")
      cat("  Features selected:", sum(glmnet::coef.glmnet(result$cvfit, s = result$lambda)[-1] != 0), 
            "out of", ncol(X), "\n")

      cat("Result for ", emb_name, " of Lasso regression: ", result$r2_test, " MSE: ", result$mse_test, "\n")
      
      results_list[[length(results_list) + 1]] <- data.frame(
        embedding_type = "baseline",
        embedding_name = emb_name,
        r2_test = result$r2_test,
        mse_test = result$mse_test,
        lambda_min = result$lambda_min,
        lambda_1se = result$lambda_1se,
        lambda_used = result$lambda,
        stringsAsFactors = FALSE
      )
      
      cat("R²=", round(result$r2_test, 4), ", MSE=", round(result$mse_test, 4), "\n", sep="")
    }
  }
  
  # Evaluate trained embeddings
  if (length(trained_emb_list) > 0) {
    cat("\n--- Trained Embeddings ---\n")
    for (emb_name in names(trained_emb_list)) {
      cat("Evaluating", emb_name, "... ")
      
      emb_data <- trained_emb_list[[emb_name]]
      # Embeddings already scaled from extraction, just remove id
      X <- as.matrix(emb_data[,c(-1)])
      

      # DEBUG: Check for zero variance columns in trained embeddings
      zero_var_cols <- apply(X, 2, function(x) {
        var_x <- var(x, na.rm=TRUE)
        is.na(var_x) || var_x == 0 || var_x < 1e-10
      })
      if (any(zero_var_cols)) {
        zero_var_indices <- which(zero_var_cols)
        cat(sprintf("[INFO] Found %d zero/low-variance columns in trained embedding %s: columns %s\n", 
                    sum(zero_var_cols),
                    emb_name,
                    paste(zero_var_indices, collapse=", ")))
        cat(sprintf("      Removing these columns before correlation computation\n"))
        X <- X[, !zero_var_cols, drop=FALSE]
      }

      # Check for missing values before running glmnet
      if (any(is.na(X)) || any(is.infinite(X))) {
        cat("SKIPPED (NA or Inf values detected in embeddings)\n")
        # ... rest of error handling
      }

      # DEBUG: Check Y variance
      Y_var <- var(Y, na.rm=TRUE)
      cat("  Y variance:", round(Y_var, 6), "\n")
      if (Y_var == 0 || is.na(Y_var)) {
        cat("  [WARNING] Y has zero variance!\n")
      }

      # DEBUG 
      cor_with_y <- apply(X, 2, function(x) {
        tryCatch({
          cor(x, Y, use="complete.obs")
        }, warning = function(w) {
          cat(sprintf("    [WARNING in cor()] Column %d: %s\n", 
                      which(apply(X, 2, identical, x))[1], 
                      conditionMessage(w)))
          return(NA)
        }, error = function(e) {
          cat(sprintf("    [ERROR in cor()] Column %d: %s\n", 
                      which(apply(X, 2, identical, x))[1], 
                      conditionMessage(e)))
          return(NA)
        })
      })
      cat("  Max |correlation| with Y:", max(abs(cor_with_y), na.rm=TRUE), "\n")
      cat("  Features with |cor| > 0.1:", sum(abs(cor_with_y) > 0.1, na.rm=TRUE), "\n")
      cat("  Features with NA correlation:", sum(is.na(cor_with_y)), "\n")
      # Check for missing values before running glmnet
      if (any(is.na(X)) || any(is.infinite(X))) {
        cat("SKIPPED (NA or Inf values detected in embeddings)\n")
        
        results_list[[length(results_list) + 1]] <- data.frame(
          embedding_type = "trained",
          embedding_name = emb_name,
          r2_test = NA,
          mse_test = NA,
          lambda_min = NA,
          lambda_1se = NA,
          lambda_used = NA,
          stringsAsFactors = FALSE
        )
        next
      }

      # Run glmnet with error handling
      result <- tryCatch({
        ridge_lasso_glmnet(
          X = X, 
          y = Y,
          alpha = alpha,
          test_prop = test_prop,
          nfolds = nfolds,
          lambda_choice = lambda_choice,
          seed = seed,
          standardize = standardize_ridge_lasso
        )
        
      }, error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        return(NULL)
      })
      
      if (is.null(result)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          embedding_type = "trained",
          embedding_name = emb_name,
          r2_test = NA,
          mse_test = NA,
          lambda_min = NA,
          lambda_1se = NA,
          lambda_used = NA,
          stringsAsFactors = FALSE
        )
        next
      }
      
      results_list[[length(results_list) + 1]] <- data.frame(
        embedding_type = "trained",
        embedding_name = emb_name,
        r2_test = result$r2_test,
        mse_test = result$mse_test,
        lambda_min = result$lambda_min,
        lambda_1se = result$lambda_1se,
        lambda_used = result$lambda,
        stringsAsFactors = FALSE
      )
      
      cat("R²=", round(result$r2_test, 4), ", MSE=", round(result$mse_test, 4), "\n", sep="")
    }
  }
  
  # ============================================================================
  # 7. Save Results
  # ============================================================================
  cat("\n=== Saving Results ===\n")
  
  # Combine results
  results_df <- do.call(rbind, results_list)
  
  # Save results CSV
  if (is_ci == "Y_age") {
    results_path <- file.path(out_dir_diagnostic, "diagnostic_results_age.csv")
  } else {
    results_path <- file.path(out_dir_diagnostic, "diagnostic_results.csv")
  }
  write.csv(results_df, results_path, row.names = FALSE)
  cat("Saved results to:", results_path, "\n")
  
  # Save diagnostic config
  diagnostic_config <- list(
    experiment_dir = experiment_dir,
    embedding_obs = embedding_obs,
    baseline_embeddings = baseline_embeddings,
    diagnostic_csv = diagnostic_csv,
    seed = seed,
    extract_trained = extract_trained,
    glmnet_params = list(
      alpha = alpha,
      test_prop = test_prop,
      nfolds = nfolds,
      lambda_choice = lambda_choice
    ),
    data_gen_params = list(
      post_non_lin = post_non_lin,
      eps_sigmaX = eps_sigmaX,
      eps_sigmaY = eps_sigmaY,
      eps_sigmaZ = eps_sigmaZ,
      embedding_orig = embedding_orig,
      confounder = confounder,
      g_z = g_z,
      is_ci = is_ci
    ),
    n_samples = length(common_ids),
    timestamp = Sys.time()
  )
  if (is_ci == "Y_age") {
    config_path <- file.path(out_dir_diagnostic, "diagnostic_config_age.json")
  } else {
    config_path <- file.path(out_dir_diagnostic, "diagnostic_config.json")
  }
  jsonlite::write_json(diagnostic_config, config_path, pretty = TRUE, auto_unbox = TRUE)
  cat("Saved config to:", config_path, "\n")
  
  cat("\n=== Diagnostic Complete ===\n")
  
  # Extract baseline results for caching
  baseline_results_only <- results_df[results_df$embedding_type == "baseline", ]
  
  # Return results
  invisible(list(
    results = results_df,
    baseline_results = baseline_results_only,
    results_path = results_path,
    config_path = config_path,
    Y_path = Y_path,
    output_dir = out_dir_diagnostic
  ))
}


# =============================================================================
# Helper Functions (extracted from pred_y.R)
# =============================================================================

#' Load confounders Z based on confounder type
#'
#' @param path_to_ukb_data Path to UKB data directory
#' @param confounder Type of confounder to load
#' @return Data frame with confounders including id column
load_Z <- function(path_to_ukb_data, confounder) {
  if (confounder == 'AS') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_Z.csv"), header=TRUE, nThread = 1)
  } else if (confounder == 'genes10') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_Z_genes10.csv"), header=TRUE, nThread = 1)
  } else if (confounder == 'ukb_z1') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z1_age.csv"), header=TRUE, nThread = 1)
  } else if (confounder == 'ukb_z2') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z2_agesex.csv"), header=TRUE, nThread = 1)
  } else if (confounder == 'ukb_z4') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z4_agesexsitesize.csv"), header=TRUE, nThread = 1)
  } else if (confounder == 'ukb_z6') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z6_agesexsitesizedateqc.csv"), header=TRUE, nThread = 1)
  } else if (confounder == 'ukb_z10') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z10_agesexsitesizedateqclocation.csv"), header=TRUE, nThread = 1)
  } else if (confounder == 'ukb_z15') {
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z15_agesexsitesizedateqcgenes.csv"), header=TRUE, nThread = 1)
  } else {
    stop("Unknown confounder type: ", confounder)
  }
  return(Z)
}

#' Check if a column is binary
#'
#' @param x Vector to check
#' @return Logical indicating if x has only 2 unique values
is_binary <- function(x) {
  unique_values <- unique(x)
  length(unique_values) == 2
}

#' Ridge/Lasso regression with glmnet
#'
#' @param X Feature matrix
#' @param y Response vector
#' @param alpha Elastic net mixing parameter (1=lasso, 0=ridge)
#' @param test_prop Proportion of data for test set
#' @param nfolds Number of CV folds
#' @param lambda_choice Lambda selection: "min" or "1se"
#' @param seed Random seed
#' @param standardize Standardize features (default: TRUE)
#' @param intercept Include intercept (default: TRUE)
#' @param return_model Return the CV model object (default: FALSE)
#'
#' @return List with test MSE, R², lambda values, and optionally the model
ridge_lasso_glmnet <- function(
  X, y,
  alpha = 1,
  test_prop = 0.2,
  nfolds = 10,
  lambda_choice = c("min", "1se"),
  seed = 42,
  standardize = TRUE,
  intercept = TRUE,
  return_model = TRUE
) {
  stopifnot(is.numeric(y))
  X <- as.matrix(X)
  stopifnot(nrow(X) == length(y))
  
  set.seed(seed)
  n <- nrow(X)
  idx_test  <- sample.int(n, size = max(1, round(test_prop * n)))
  idx_train <- setdiff(seq_len(n), idx_test)
  
  Xtr <- X[idx_train, , drop = FALSE]
  ytr <- y[idx_train]
  Xte <- X[idx_test,  , drop = FALSE]
  yte <- y[idx_test]
  
  # Cross-validated glmnet
  suppressPackageStartupMessages(require(glmnet))
  cvfit <- glmnet::cv.glmnet(
    x = Xtr, y = ytr,
    family = "gaussian",
    alpha = alpha,
    nfolds = nfolds,
    standardize = standardize,
    intercept = intercept
  )

  # After cvfit is created, add:
  cat("  CV error curve:\n")
  cat("    Lambda range: [", min(cvfit$lambda), ", ", max(cvfit$lambda), "]\n")
  cat("    CV error range: [", min(cvfit$cvm), ", ", max(cvfit$cvm), "]\n")
  cat("    Lambda.min index:", which(cvfit$lambda == cvfit$lambda.min), 
    "out of", length(cvfit$lambda), "\n")
  
  lambda_choice <- match.arg(lambda_choice)
  s_use <- if (lambda_choice == "min") cvfit$lambda.min else cvfit$lambda.1se
  
  pred <- drop(predict(cvfit, newx = Xte, s = s_use))
  mse  <- mean((yte - pred)^2)
  r2   <- 1 - sum((yte - pred)^2) / sum((yte - mean(yte))^2)
  
  out <- list(
    mse_test = as.numeric(mse),
    r2_test  = as.numeric(r2),
    lambda   = s_use,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    alpha = alpha,
    idx_train = idx_train,
    idx_test = idx_test,
    y_test = yte,
    yhat_test = pred
  )
  if (return_model) out$cvfit <- cvfit
  out
}
