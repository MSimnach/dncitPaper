### After fitting networks, extracting embeddings, running diagnostics
# 1. Aim: compare if R^2 of embeddings is larger for large FAST ROIs
# Select ROIs with large regions and check if for simulation runs their weights were large
# Expectation: for these runs, R^2 should be larger
# 2. Aim: check if trained embeddings profit more from large FAST ROIs compared to baseline embeddings
devtools::load_all()

ids_fast <- data.table::fread("/sc/home/marco.simnacher/dncitPaper/inst/extdata/ids/ids_fastsurfer.csv", header=TRUE, nThread = 1)
ukb_fast <- data.table::fread("/sc/home/marco.simnacher/ukbiobank/data/ukb_fastsurfer.csv", header=TRUE, nThread = 1)
#substring of colnames of ukb_fast before -2.0
colnames_fast <- substring(colnames(ukb_fast)[-c(1)], 1, nchar(colnames(ukb_fast)[-c(1)]) - 4)
all(colnames_fast== as.character(ids_fast$ID))

large_fast_rois_idx <- c(
  5, 6, 7,                          # Brain-stem, Caudate L/R
  14, 15, 16, 17, 18, 19,           # Crus I & II (L/R/vermis)
  32, 33,                           # Hippocampus L/R
  34, 35,                           # I–IV Cerebellum L/R
  55, 56, 57, 58,                   # Lateral Occipital (inf/sup, L/R)
  59, 60,                           # Lingual Gyrus L/R
  61, 62,                           # Middle Frontal Gyrus L/R
  63, 64, 65, 66, 67, 68,           # Middle Temporal (ant/post/TO, L/R)
  69, 70,                           # Occipital Fusiform Gyrus L/R
  71, 72,                           # Occipital Pole L/R
  73, 74,                           # Pallidum L/R
  87, 88,                           # Postcentral Gyrus L/R
  89, 90,                           # Precentral Gyrus L/R
  93, 94,                           # Putamen L/R
  97, 98,                           # Superior Frontal Gyrus L/R
  99, 100,                          # Superior Parietal Lobule L/R
  101, 102, 103, 104,               # Superior Temporal (ant/post, L/R)
  107, 108, 109, 110,               # Supramarginal (ant/post, L/R)
  111, 112, 113, 114,               # Temporal Fusiform (ant/post, L/R)
  119, 120,                         # Thalamus L/R
  121, 122,                         # V Cerebellum L/R
  123, 124, 125,                    # VI Cerebellum (L/R/vermis)
  126, 127, 128,                    # VIIIa Cerebellum (L/R/vermis)
  129, 130, 131,                    # VIIIb Cerebellum (L/R/vermis)
  132, 133, 134,                    # VIIb Cerebellum (L/R/vermis)
  135, 136                          # Ventral Striatum L/R
)
# smaller ROI set
large_fast_rois_idx <- c(
  # Subcortical
  5, 6,      # Caudate L/R
  93, 94,    # Putamen L/R
  73, 74,    # Pallidum L/R
  119, 120,  # Thalamus L/R
  7,         # Brainstem
  
  # Motor / Parietal cortex
  89, 90,    # Precentral
  87, 88,    # Postcentral
  99, 100,   # Superior Parietal Lobule
  
  # Occipital regions
  55, 56, 57, 58,  # Lat Occipital
  59, 60,          # Lingual
  71, 72,          # Occipital Pole
  
  # Cerebellum
  121, 122,        # V
  123, 124, 125    # VI (L/R+vermis)
)
large_fast_rois_ids <- colnames_fast[large_fast_rois_idx]
ids_fast[as.character(ids_fast$ID) %in% large_fast_rois_ids,]

# for all seeds, get weights for X_orig in data gen
seeds <- c(601:604)
debug_Y <- TRUE
eps_sigmaY <- 0.5
post_non_lin <- 1
g_z <- "squared"
xz_mode <- "Sigma=I_p"
Z <- matrix(rnorm(100 * 2), 100, 2)
X_orig <- matrix(rnorm(100 * 139), 100, 139)
#save large share in df 
seeds_large_rois <- data.frame(seed = seeds, share = NA, n_large_active = NA)
weights_X_orig_list <- list()
weights_X_orig_large_list <- list()
# compute also ratio of largest weights compared to all other weights
weights_X_orig_large_ratio_list <- list()
for(seed in seeds){
    set.seed(seed + 999999)
  Y_info <- y_from_xz(Z, eps_sigmaY, X=X_orig, post_non_lin=post_non_lin, g_z=g_z, xz_mode=xz_mode, debug=debug_Y)
  weights_X_orig <- Y_info$info$beta_X
  wx_large <- weights_X_orig[large_fast_rois_idx]
  weights_X_orig_list[[as.character(seed)]] <- weights_X_orig
  weights_X_orig_large_list[[as.character(seed)]] <- weights_X_orig[large_fast_rois_idx]
  seeds_large_rois$share[seeds_large_rois$seed == seed] <- sum(abs(wx_large)) / sum(abs(weights_X_orig))
  seeds_large_rois$n_large_active[seeds_large_rois$seed == seed] <- sum(abs(wx_large) > 0)
  weights_X_orig_large_ratio_list[[as.character(seed)]] <- max(abs(weights_X_orig)) / sum(abs(weights_X_orig))
}

seeds_large_rois$share_rank <- rank(seeds_large_rois$share)
seeds_large_rois$weights_ratio <- unlist(weights_X_orig_large_ratio_list)
seeds_large_rois$weights_ratio_rank <- rank(unlist(seeds_large_rois$weights_ratio))

cat("=== Collecting Diagnostic Results ===\n")

# Parameters
n_samples <- c(460, 1100, 5000, 10000)#, 10000)#c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)
eps_sigmaY_all <- c(0.5)
conditions <- c("CI", "No_CI")
base_paths <- list(
  "CI" = "/sc/home/marco.simnacher/ukbiobank/data/CI",
  "No_CI" = "/sc/home/marco.simnacher/ukbiobank/data/No_CI"
)

# Initialize lists to store results for Y and age
all_results_y <- list()
all_results_age <- list()

# Loop through all combinations and read results
for (condition in conditions) {
  base_path <- base_paths[[condition]]
  for (n_sample in n_samples) {
    for (seed in seeds) {
      for (eps_sigmaY in eps_sigmaY_all) {
        diagnostic_dir <- file.path(
          base_path, 
          n_sample, 
          seed, 
          paste0("eps_sigmaY=", eps_sigmaY),
          "diagnostic"
        )
        
        # Read Y diagnostic results
        results_path_y <- file.path(diagnostic_dir, "diagnostic_results.csv")
        if (file.exists(results_path_y)) {
          results_df <- fread(results_path_y, nThread = 1)
          results_df$n_sample <- n_sample
          results_df$seed <- seed
          results_df$eps_sigmaY <- eps_sigmaY
          results_df$condition <- condition
          all_results_y[[length(all_results_y) + 1]] <- results_df
        } else {
          warning("File not found: ", results_path_y)
        }
        
        # Read age diagnostic results
        results_path_age <- file.path(diagnostic_dir, "diagnostic_results_age.csv")
        if (file.exists(results_path_age)) {
          results_df <- fread(results_path_age, nThread = 1)
          results_df$n_sample <- n_sample
          results_df$seed <- seed
          results_df$eps_sigmaY <- eps_sigmaY
          results_df$condition <- condition
          all_results_age[[length(all_results_age) + 1]] <- results_df
        } else {
          warning("File not found: ", results_path_age)
        }
      }
    }
  }
}


combined_results_y <- rbindlist(all_results_y)
constant_embeddings <- c("fastsurfer", "freesurfer", "condVAE", "medicalnet")#, "pooled_brainsynth", "tucker_brainsynth")
varying_embeddings <- c('medicalnet_ft', 'scratch')
embeddings <- c(constant_embeddings, varying_embeddings)


# correlation R^2 per embedding and large share rank
cor_results <- data.frame(embedding = embeddings, r2_share = NA, test_p_val = NA, r2_share_rank = NA, r2_scratch_medicalnet = NA)
# Create a list of matrices for R^2 values, one matrix per sample size and condition
# Rows: seeds, Columns: embeddings
r2_matrices_by_sample_size_and_condition <- list()

for (cond in conditions) {  # Changed from 'condition' to 'cond'
  r2_matrices_by_sample_size_and_condition[[cond]] <- list()
  
  for (n_samp in n_samples) {  # Changed from 'n_sample' to 'n_samp'
    # Use explicit variable names that don't match column names
    data_subset <- combined_results_y[n_sample == n_samp & condition == cond]
    
    # Debug: check how many rows we got
    cat(sprintf("Condition: %s, n_sample: %d, rows: %d\n", cond, n_samp, nrow(data_subset)))
    
    r2_matrix <- matrix(NA, 
                        nrow = length(seeds), 
                        ncol = length(embeddings),
                        dimnames = list(as.character(seeds), embeddings))
    
    for (emb in embeddings) {  # Changed from 'embedding' to 'emb'
      emb_data <- data_subset[embedding_name == emb]
      # Match seeds to fill in correct rows
      matched_idx <- match(emb_data$seed, seeds)
      r2_matrix[matched_idx, emb] <- emb_data$r2_test
    }
    
    r2_matrices_by_sample_size_and_condition[[cond]][[as.character(n_samp)]] <- r2_matrix
  }
}

# View structure
names(r2_matrices_by_sample_size_and_condition)
# Example: view first few rows and columns of the first matrix
r2_matrices_by_sample_size_and_condition[["No_CI"]][["10000"]] 
r2_matrices_by_sample_size_and_condition[["CI"]][["460"]] 

# Compute correlation matrices between embeddings for each sample size and condition
embedding_correlations <- list()

for (condition in conditions) {
  embedding_correlations[[condition]] <- list()
  
  for (n_sample in n_samples) {
    r2_matrix <- r2_matrices_by_sample_size_and_condition[[condition]][[as.character(n_sample)]]
    
    # Compute correlation between columns (embeddings)
    # use = "pairwise.complete.obs" handles any missing values
    cor_matrix <- cor(r2_matrix, method = "pearson", use = "pairwise.complete.obs")
    
    embedding_correlations[[condition]][[as.character(n_sample)]] <- cor_matrix
  }
}

# View results
# Example: correlation between embeddings for n=10000, No_CI condition
print("Correlations for n=10000, No_CI:")
print(round(embedding_correlations[["No_CI"]][["10000"]], 3))

# Example: correlation between embeddings for n=460, CI condition
print("\nCorrelations for n=460, CI:")
print(round(embedding_correlations[["CI"]][["460"]], 3))

# If you want to see how correlations change with sample size for a specific pair:
# Extract correlation between two specific embeddings across sample sizes
print("\nCorrelation between fastsurfer and scratch across sample sizes (No_CI):")
for (n_sample in n_samples) {
  cor_val <- embedding_correlations[["No_CI"]][[as.character(n_sample)]]["fastsurfer", "scratch"]
  cat(sprintf("n=%5d: %.3f\n", n_sample, cor_val))
}


# Compute comprehensive correlation statistics for each sample size, condition, and embedding
# Create matrices with embeddings as rows and different correlation metrics as columns
cor_results_matrices <- list()

for (cond in conditions) {
  cor_results_matrices[[cond]] <- list()
  
  for (n_samp in n_samples) {
    # Initialize matrix: rows = embeddings, columns = different correlation metrics
    cor_matrix <- matrix(NA, 
                         nrow = length(embeddings), 
                         ncol = 4,
                         dimnames = list(embeddings, 
                                       c("r2_share", "r2_share_rank", "p_value", "r2_weights_ratio")))
    
    # Get R^2 matrix for this sample size and condition
    r2_matrix <- r2_matrices_by_sample_size_and_condition[[cond]][[as.character(n_samp)]]
    
    for (emb in embeddings) {
      # Get R^2 values for this embedding across seeds
      r2_vals <- r2_matrix[, emb]
      
      # 1. Correlation with share of large ROIs
      cor_matrix[emb, "r2_share"] <- cor(r2_vals, seeds_large_rois$share, use = "complete.obs")
      
      # 2. Correlation with rank of share
      cor_matrix[emb, "r2_share_rank"] <- cor(r2_vals, seeds_large_rois$share_rank, use = "complete.obs")
      
      # 3. P-value from cor.test (testing correlation with share)
      cor_test_result <- cor.test(r2_vals, seeds_large_rois$share)
      cor_matrix[emb, "p_value"] <- cor_test_result$p.value
      
      # 4. Correlation with weights_ratio
      cor_matrix[emb, "r2_weights_ratio"] <- cor(r2_vals, unlist(seeds_large_rois$weights_ratio), use = "complete.obs")
    }
    
    cor_results_matrices[[cond]][[as.character(n_samp)]] <- cor_matrix
  }
}

# Print results for each condition and sample size
cat("\n=== Correlation Results: R^2 vs Large ROI Metrics ===\n")
for (cond in conditions) {
  cat(sprintf("\n\n### Condition: %s ###\n", cond))
  for (n_samp in n_samples) {
    cat(sprintf("\nn_sample = %d:\n", n_samp))
    print(round(cor_results_matrices[[cond]][[as.character(n_samp)]], 4))
  }
}

# Create a combined data frame for easier analysis and plotting
cor_results_df_list <- list()
for (cond in conditions) {
  for (n_samp in n_samples) {
    cor_mat <- cor_results_matrices[[cond]][[as.character(n_samp)]]
    df <- as.data.frame(cor_mat)
    df$embedding <- rownames(df)
    df$condition <- cond
    df$n_sample <- n_samp
    cor_results_df_list[[length(cor_results_df_list) + 1]] <- df
  }
}
cor_results_df <- rbindlist(cor_results_df_list)

# Reorder columns for readability
cor_results_df <- cor_results_df[, c("condition", "n_sample", "embedding", "r2_share", "r2_share_rank", "p_value", "r2_weights_ratio")]

cat("\n\n=== Combined Correlation Results (Data Frame) ===\n")
print(cor_results_df)

# Summary: Compare trained vs baseline embeddings
cat("\n\n=== Summary: Trained vs Baseline Embeddings ===\n")
cat("\nMean correlation (r2_share) across sample sizes:\n")
summary_by_embedding <- cor_results_df[, .(
  mean_r2_share = mean(r2_share, na.rm = TRUE),
  mean_r2_share_rank = mean(r2_share_rank, na.rm = TRUE),
  mean_r2_weights_ratio = mean(r2_weights_ratio, na.rm = TRUE),
  min_p_value = min(p_value, na.rm = TRUE),
  median_p_value = median(p_value, na.rm = TRUE)
), by = .(embedding, condition)]
print(summary_by_embedding[order(condition, -mean_r2_share)])

# Check if correlations improve with sample size
cat("\n\n=== Correlation Trends by Sample Size (No_CI condition) ===\n")
trend_analysis <- cor_results_df[condition == "No_CI", .(embedding, n_sample, r2_share, p_value)]
trend_analysis_wide <- dcast(trend_analysis, embedding ~ n_sample, value.var = "r2_share")
print(trend_analysis_wide)

# Identify embeddings with significant positive correlations
cat("\n\n=== Embeddings with Significant Positive Correlation (p < 0.05) ===\n")
significant_results <- cor_results_df[p_value < 0.05 & r2_share > 0]
print(significant_results[order(condition, n_sample, -r2_share)])

### Same for p-values of CITs and large shares
cat("\n\n=== Loading P-values from CIT Tests ===\n")

# Define p-values directories
pval_base_paths <- list(
  "CI" = "/sc/home/marco.simnacher/dncitPaper/Results/CI/p-values/seeds_601_625",
  "No_CI" = "/sc/home/marco.simnacher/dncitPaper/Results/No_CI/p-values/seeds_601_625"
)

# Function to extract embedding name from filename
extract_embedding_from_filename <- function(filename) {
  # Pattern: 1_0_[eps]_0_fastsurfer_[embedding]_ukb_z1_squared_RCOT_1.csv
  # Extract everything between "fastsurfer_" and "_ukb"
  pattern <- "fastsurfer_(.+?)_ukb"
  match <- regmatches(filename, regexec(pattern, filename))
  if (length(match[[1]]) > 1) {
    return(match[[1]][2])
  } else {
    return(NA)
  }
}

# Load all p-value CSV files
pval_data_list <- list()

for (cond in conditions) {
  pval_dir <- pval_base_paths[[cond]]
  
  # Get all CSV files in the directory
  csv_files <- list.files(pval_dir, pattern = "\\.csv$", full.names = TRUE)
  
  cat(sprintf("\nCondition: %s - Found %d CSV files\n", cond, length(csv_files)))
  
  for (csv_file in csv_files) {
    filename <- basename(csv_file)
    embedding_name <- extract_embedding_from_filename(filename)
    
    if (is.na(embedding_name)) {
      warning(sprintf("Could not extract embedding name from: %s", filename))
      next
    }
    
    cat(sprintf("  Loading: %s (embedding: %s)\n", filename, embedding_name))
    
    # Read the CSV
    pval_df <- fread(csv_file, nThread = 1)
    
    # Drop the first column by position (contains just 1:25 row indices)
    # fread assigns it "V1" which conflicts with actual data column names
    pval_df <- pval_df[, -1]
    
    # Now we have columns for the 4 sample sizes (460, 1100, 5000, 10000)
    # Rename columns to avoid confusion
    sample_size_cols <- names(pval_df)[1:4]
    
    # Add seed information based on row number
    pval_df$seed_idx <- 1:nrow(pval_df)
    pval_df$seed <- seeds[pval_df$seed_idx]
    
    # Reshape to long format using the actual column names
    pval_long <- melt(pval_df, 
                      id.vars = c("seed", "seed_idx"),
                      measure.vars = sample_size_cols,
                      variable.name = "sample_size_col",
                      value.name = "p_value")
    
    # Map V1-V4 to actual sample sizes
    pval_long$n_sample <- n_samples[as.integer(gsub("V", "", pval_long$sample_size_col))]
    
    # Add metadata
    pval_long$embedding <- embedding_name
    pval_long$condition <- cond
    
    # Store in list
    pval_data_list[[length(pval_data_list) + 1]] <- pval_long[, .(seed, n_sample, embedding, condition, p_value)]
  }
}

# Combine all p-value data
combined_pvals <- rbindlist(pval_data_list)

cat(sprintf("\nTotal p-value records loaded: %d\n", nrow(combined_pvals)))
cat(sprintf("Unique embeddings: %d\n", length(unique(combined_pvals$embedding))))
cat(sprintf("Embeddings found: %s\n", paste(unique(combined_pvals$embedding), collapse = ", ")))

# Create matrices for p-values similar to R^2 matrices
# Rows: seeds, Columns: embeddings
pval_matrices_by_sample_size_and_condition <- list()
available_embeddings <- unique(combined_pvals$embedding)

for (cond in conditions) {
  pval_matrices_by_sample_size_and_condition[[cond]] <- list()
  
  for (n_samp in n_samples) {
    # Filter data for this condition and sample size
    data_subset <- combined_pvals[n_sample == n_samp & condition == cond]
    
    # Debug info
    cat(sprintf("P-values - Condition: %s, n_sample: %d, rows: %d\n", cond, n_samp, nrow(data_subset)))
    
    # Create matrix
    pval_matrix <- matrix(NA, 
                          nrow = length(seeds), 
                          ncol = length(available_embeddings),
                          dimnames = list(as.character(seeds), available_embeddings))
    
    for (emb in available_embeddings) {
      emb_data <- data_subset[embedding == emb]
      if (nrow(emb_data) > 0) {
        # Match seeds to fill in correct rows
        matched_idx <- match(emb_data$seed, seeds)
        pval_matrix[matched_idx, emb] <- emb_data$p_value
      }
    }
    
    pval_matrices_by_sample_size_and_condition[[cond]][[as.character(n_samp)]] <- pval_matrix
  }
}

# Display example matrices
cat("\n=== Example P-value Matrices ===\n")
if (!is.null(pval_matrices_by_sample_size_and_condition[["No_CI"]][["10000"]])) {
  cat("\nNo_CI, n=10000 (first 5 seeds):\n")
  print(pval_matrices_by_sample_size_and_condition[["No_CI"]][["10000"]][1:5, , drop = FALSE])
}
if (!is.null(pval_matrices_by_sample_size_and_condition[["CI"]][["460"]])) {
  cat("\nCI, n=460 (first 5 seeds):\n")
  print(pval_matrices_by_sample_size_and_condition[["CI"]][["460"]][1:5, , drop = FALSE])
}

# Compute correlations between p-values and large ROI metrics
pval_cor_results_matrices <- list()

for (cond in conditions) {
  pval_cor_results_matrices[[cond]] <- list()
  
  for (n_samp in n_samples) {
    # Get p-value matrix for this sample size and condition
    pval_matrix <- pval_matrices_by_sample_size_and_condition[[cond]][[as.character(n_samp)]]
    
    if (is.null(pval_matrix) || all(is.na(pval_matrix))) {
      cat(sprintf("Warning: No p-value data for condition=%s, n_sample=%d\n", cond, n_samp))
      next
    }
    
    # Initialize correlation matrix
    pval_cor_matrix <- matrix(NA, 
                              nrow = length(available_embeddings), 
                              ncol = 4,
                              dimnames = list(available_embeddings, 
                                            c("pval_share", "pval_share_rank", "p_value_cortest", "pval_weights_ratio")))
    
    for (emb in available_embeddings) {
      # Get p-values for this embedding across seeds
      pvals <- pval_matrix[, emb]
      
      # Skip if all NA
      if (all(is.na(pvals))) next
      
      # 1. Correlation with share of large ROIs
      pval_cor_matrix[emb, "pval_share"] <- cor(pvals, seeds_large_rois$share, use = "complete.obs")
      
      # 2. Correlation with rank of share
      pval_cor_matrix[emb, "pval_share_rank"] <- cor(pvals, seeds_large_rois$share_rank, use = "complete.obs")
      
      # 3. P-value from cor.test (testing correlation with share)
      cor_test_result <- cor.test(pvals, seeds_large_rois$share)
      pval_cor_matrix[emb, "p_value_cortest"] <- cor_test_result$p.value
      
      # 4. Correlation with weights_ratio
      pval_cor_matrix[emb, "pval_weights_ratio"] <- cor(pvals, unlist(seeds_large_rois$weights_ratio), use = "complete.obs")
    }
    
    pval_cor_results_matrices[[cond]][[as.character(n_samp)]] <- pval_cor_matrix
  }
}

# Print correlation results for p-values
cat("\n\n=== Correlation Results: P-values vs Large ROI Metrics ===\n")
for (cond in conditions) {
  cat(sprintf("\n\n### Condition: %s ###\n", cond))
  for (n_samp in n_samples) {
    pval_cor_mat <- pval_cor_results_matrices[[cond]][[as.character(n_samp)]]
    if (!is.null(pval_cor_mat)) {
      cat(sprintf("\nn_sample = %d:\n", n_samp))
      print(round(pval_cor_mat, 4))
    }
  }
}

# Create a combined data frame for p-value correlations
pval_cor_results_df_list <- list()
for (cond in conditions) {
  for (n_samp in n_samples) {
    pval_cor_mat <- pval_cor_results_matrices[[cond]][[as.character(n_samp)]]
    if (!is.null(pval_cor_mat)) {
      df <- as.data.frame(pval_cor_mat)
      df$embedding <- rownames(df)
      df$condition <- cond
      df$n_sample <- n_samp
      pval_cor_results_df_list[[length(pval_cor_results_df_list) + 1]] <- df
    }
  }
}

if (length(pval_cor_results_df_list) > 0) {
  pval_cor_results_df <- rbindlist(pval_cor_results_df_list, fill = TRUE)
  
  # Reorder columns for readability
  pval_cor_results_df <- pval_cor_results_df[, c("condition", "n_sample", "embedding", "pval_share", "pval_share_rank", "p_value_cortest", "pval_weights_ratio")]
  
  cat("\n\n=== Combined P-value Correlation Results (Data Frame) ===\n")
  print(pval_cor_results_df)
  
  # Summary: Compare embeddings
  cat("\n\n=== Summary: P-value Correlations by Embedding ===\n")
  cat("\nMean correlation (pval_share) across sample sizes:\n")
  pval_summary_by_embedding <- pval_cor_results_df[, .(
    mean_pval_share = mean(pval_share, na.rm = TRUE),
    mean_pval_share_rank = mean(pval_share_rank, na.rm = TRUE),
    mean_pval_weights_ratio = mean(pval_weights_ratio, na.rm = TRUE),
    min_p_value_cortest = min(p_value_cortest, na.rm = TRUE),
    median_p_value_cortest = median(p_value_cortest, na.rm = TRUE)
  ), by = .(embedding, condition)]
  print(pval_summary_by_embedding[order(condition, -mean_pval_share)])
  
  # Identify embeddings with significant negative correlations (lower p-values with larger ROI share)
  cat("\n\n=== Embeddings with Significant Negative Correlation (p < 0.05, expecting lower p-values) ===\n")
  significant_pval_results <- pval_cor_results_df[p_value_cortest < 0.05 & pval_share < 0]
  if (nrow(significant_pval_results) > 0) {
    print(significant_pval_results[order(condition, n_sample, pval_share)])
  } else {
    cat("No significant negative correlations found.\n")
  }
  
  # Check if correlations change with sample size
  cat("\n\n=== P-value Correlation Trends by Sample Size ===\n")
  for (cond in conditions) {
    cat(sprintf("\nCondition: %s\n", cond))
    pval_trend_analysis <- pval_cor_results_df[condition == cond, .(embedding, n_sample, pval_share, p_value_cortest)]
    pval_trend_analysis_wide <- dcast(pval_trend_analysis, embedding ~ n_sample, value.var = "pval_share")
    print(pval_trend_analysis_wide)
  }
} else {
  cat("\n\nNo p-value correlation results available.\n")
}

# Compute correlation between R^2 and p-values for each embedding, sample size, and condition
cat("\n\n=== Correlation between R^2 and P-values ===\n")

r2_pval_correlations <- list()

for (cond in conditions) {
  if (cond == "CI") {
    next
  }
  for (n_samp in n_samples) {
    # Get R^2 and p-value matrices
    r2_mat <- r2_matrices_by_sample_size_and_condition[[cond]][[as.character(n_samp)]]
    pval_mat <- pval_matrices_by_sample_size_and_condition[[cond]][[as.character(n_samp)]]
    
    if (is.null(r2_mat) || is.null(pval_mat)) next
    
    # Find common embeddings between both matrices
    common_embeddings <- intersect(colnames(r2_mat), colnames(pval_mat))
    
    for (emb in common_embeddings) {
      r2_vals <- r2_mat[, emb]
      pval_vals <- pval_mat[, emb]
      
      # Compute correlation (we expect negative: higher R^2 -> lower p-value)
      cor_result <- cor.test(r2_vals, pval_vals, use = "complete.obs")
      
      r2_pval_correlations[[length(r2_pval_correlations) + 1]] <- data.frame(
        condition = cond,
        n_sample = n_samp,
        embedding = emb,
        correlation = cor_result$estimate,
        p_value = cor_result$p.value,
        n_obs = sum(!is.na(r2_vals) & !is.na(pval_vals))
      )
    }
  }
}

# Combine into data frame
r2_pval_cor_df <- rbindlist(r2_pval_correlations)
rownames(r2_pval_cor_df) <- NULL

cat("\n=== R^2 vs P-value Correlations ===\n")
print(r2_pval_cor_df[order(embedding, condition, n_sample)])

cat("\n\n=== Summary: Mean correlation by embedding ===\n")
r2_pval_summary <- r2_pval_cor_df[, .(
  mean_correlation = mean(correlation, na.rm = TRUE),
  median_correlation = median(correlation, na.rm = TRUE),
  min_p_value = min(p_value, na.rm = TRUE)
), by = .(embedding, condition)]
print(r2_pval_summary[order(condition, mean_correlation)])

cat("\n\n=== Significant negative correlations (higher R^2 -> lower p-value) ===\n")
significant_r2_pval <- r2_pval_cor_df[p_value < 0.05 & correlation < 0]
if (nrow(significant_r2_pval) > 0) {
  print(significant_r2_pval[order(condition, n_sample, correlation)])
} else {
  cat("No significant negative correlations found.\n")
}

# ============================================================================
# Create Boxplots of P-values
# ============================================================================
cat("\n\n=== Creating P-value Boxplots ===\n")

library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# Prepare data for plotting
if (nrow(combined_pvals) > 0) {
  # Define embedding types
  constant_embeddings_plot <- intersect(constant_embeddings, available_embeddings)
  varying_embeddings_plot <- intersect(varying_embeddings, available_embeddings)
  
  cat("Constant embeddings available:", paste(constant_embeddings_plot, collapse = ", "), "\n")
  cat("Varying embeddings available:", paste(varying_embeddings_plot, collapse = ", "), "\n")
  
  # Set minimum p-value for plotting (to avoid log10(0) = -Inf)
  min_pval_plot <- 1e-16
  
  # Prepare plot data and replace 0 or extremely small p-values
  plot_data_pval <- combined_pvals %>%
    mutate(
      p_value = ifelse(p_value < min_pval_plot | p_value == 0, min_pval_plot, p_value),
      n_sample = factor(n_sample, levels = n_samples),
      condition = factor(condition, levels = conditions),
      embedding = factor(embedding)
    )
  
  # Report how many p-values were adjusted
  n_adjusted <- sum(combined_pvals$p_value < min_pval_plot | combined_pvals$p_value == 0, na.rm = TRUE)
  if (n_adjusted > 0) {
    cat(sprintf("Note: %d p-values (%.1f%%) were < %g or zero, set to %g for plotting\n", 
                n_adjusted, 100*n_adjusted/nrow(combined_pvals), min_pval_plot, min_pval_plot))
  }
  
  # Create color palette for all embeddings
  all_embeddings_plot <- c(constant_embeddings_plot, varying_embeddings_plot)
  n_colors <- length(all_embeddings_plot)
  base_colors <- RColorBrewer::brewer.pal(min(8, n_colors), "Set2")
  color_palette_pval <- setNames(
    colorRampPalette(base_colors)(n_colors),
    all_embeddings_plot
  )
  
  # Calculate y-axis limits (log scale) - use adjusted data so no Inf values
  y_limits_pval <- c(min_pval_plot, 1)
  
  # Combine constant and varying embeddings for plotting across all sample sizes
  plot_data_pval_all <- plot_data_pval %>%
    filter(embedding %in% all_embeddings_plot, condition == "No_CI") %>%
    mutate(embedding = factor(embedding, levels = all_embeddings_plot))
  
  # Single panel: All embeddings across sample sizes
  p_pval_combined <- plot_data_pval_all %>%
    ggplot(aes(x = n_sample, y = p_value, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
    labs(
      x = "Sample Size",
      y = "P-value",
      fill = "Embedding",
      title = expression("CIT P-values: No_CI Condition (" * epsilon[sigma[Y]] * " = 0.5)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    scale_fill_manual(values = color_palette_pval) +
    scale_y_log10(
      breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.5, 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = y_limits_pval
    )
  
  # Save plot
  output_dir <- "/sc/home/marco.simnacher/dncitPaper/inst/diagnostic_embedding/figures"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\nSaving p-value boxplots:\n")
  
  # Save as PNG
  png_path_pval <- file.path(output_dir, paste0("pvalues_boxplot_", min(seeds), "_", max(seeds), ".png"))
  ggsave(
    png_path_pval, 
    plot = p_pval_combined, 
    width = 14, 
    height = 7, 
    units = "in", 
    dpi = 300
  )
  cat("  Saved PNG:", png_path_pval, "\n")
  
  # Save as PDF
  pdf_path_pval <- file.path(output_dir, paste0("pvalues_boxplot_", min(seeds), "_", max(seeds), ".pdf"))
  ggsave(
    pdf_path_pval, 
    plot = p_pval_combined, 
    width = 14, 
    height = 7, 
    units = "in"
  )
  cat("  Saved PDF:", pdf_path_pval, "\n")
  
  cat("\n=== P-value Visualization Complete ===\n")
} else {
  cat("No p-value data available for plotting.\n")
}

