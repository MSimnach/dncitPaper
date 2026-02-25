### Correlation between CIT P-values and Diagnostics
# Compute correlations between CIT p-values and diagnostic metrics
# (R^2 original, R^2 residual, F-test, globaltest, dcor.test p-values)

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(paletteer)

# Configuration
seeds <- c(1:100)
n_samples <- c(256, 460, 825, 1100, 5000)#, 10000)
conditions <- c("CI", "No_CI")
eps_sigmaY <- 0.5

# Helper functions
extract_embedding_from_filename <- function(filename) {
  pattern <- "fastsurfer_(.+?)_ukb"
  match <- regmatches(filename, regexec(pattern, filename))
  if (length(match[[1]]) > 1) return(match[[1]][2]) else return(NA)
}

extract_cit_from_filename <- function(filename) {
  pattern <- "squared_(.+?)\\.csv"
  match <- regmatches(filename, regexec(pattern, filename))
  if (length(match[[1]]) > 1) return(match[[1]][2]) else return(NA)
}

# ============================================================================
# Load CIT P-values
# ============================================================================
cat("=== Loading CIT P-values ===\n")

results_base_path <- "/sc/home/marco.simnacher/dncitPaper/Results"
pval_data_list <- list()

for (cond in conditions) {
  pval_dir <- file.path(results_base_path, cond, "p-values", sprintf("seeds_%d_%d", min(seeds), max(seeds)))
  
  if (!dir.exists(pval_dir)) {
    warning(sprintf("Directory not found: %s", pval_dir))
    next
  }
  
  # Filter CSV files by eps_sigmaY pattern (e.g., 1_0_0.5_0_...)
  pattern_with_eps <- sprintf("^1_0_%g_0_.*\\.csv$", eps_sigmaY)
  csv_files <- list.files(pval_dir, pattern = pattern_with_eps, full.names = TRUE)
  cat(sprintf("Condition: %s - Found %d CSV files for eps_sigmaY=%s\n", cond, length(csv_files), eps_sigmaY))
  
  for (csv_file in csv_files) {
    filename <- basename(csv_file)
    embedding_name <- extract_embedding_from_filename(filename)
    cit_name <- extract_cit_from_filename(filename)
    
    if (is.na(embedding_name) || is.na(cit_name)) next
    
    pval_df <- fread(csv_file, nThread = 1)
    pval_df <- pval_df[, -1]  # Remove row index column
    
    sample_size_cols <- names(pval_df)[c(2,4,5)] #adapt to sample sizes
    pval_df$seed_idx <- 1:nrow(pval_df)
    pval_df$seed <- seeds[pval_df$seed_idx]
    
    pval_long <- melt(pval_df, 
                      id.vars = c("seed", "seed_idx"),
                      measure.vars = sample_size_cols,
                      variable.name = "sample_size_col",
                      value.name = "p_value")
    
    pval_long$n_sample <- n_samples[as.integer(gsub("V", "", pval_long$sample_size_col))]
    pval_long$embedding <- embedding_name
    pval_long$condition <- cond
    pval_long$cit <- cit_name
    
    pval_data_list[[length(pval_data_list) + 1]] <- pval_long[, .(seed, n_sample, embedding, condition, cit, p_value)]
  }
}

combined_pvals <- rbindlist(pval_data_list)
cat(sprintf("Total CIT p-value records: %d\n", nrow(combined_pvals)))
cat(sprintf("Unique CITs: %s\n", paste(unique(combined_pvals$cit), collapse = ", ")))
cat(sprintf("Unique embeddings: %s\n\n", paste(unique(combined_pvals$embedding), collapse = ", ")))

n_samples <- c(460,1100,5000)
# ============================================================================
# Load Diagnostic Results
# ============================================================================
cat("=== Loading Diagnostic Results ===\n")

base_paths <- list(
  "CI" = "/sc/home/marco.simnacher/ukbiobank/data/CI",
  "No_CI" = "/sc/home/marco.simnacher/ukbiobank/data/No_CI"
)

all_diagnostics <- list()

for (condition in conditions) {
  base_path <- base_paths[[condition]]
  for (n_sample in n_samples) {
    for (seed in seeds) {
      diagnostic_path <- file.path(
        base_path, n_sample, seed, 
        paste0("eps_sigmaY=", eps_sigmaY),
        "diagnostic/diagnostic_results.csv"
      )
      
      if (file.exists(diagnostic_path)) {
        diag_df <- fread(diagnostic_path, nThread = 1)
        
        # Check if all required columns exist
        required_cols <- c("embedding_name", "y_type", "r2_test", 
                          "f_test_pvalue", "globaltest_pvalue", "dcor_pvalue",
                          "pcm_pvalue", "rcot_pvalue", "pcm_pvalue_split", "rcot_pvalue_split")
        if (!all(required_cols %in% names(diag_df))) {
          warning(sprintf("Skipping %s - missing required columns", diagnostic_path))
          next
        }
        
        diag_df$n_sample <- n_sample
        diag_df$seed <- seed
        diag_df$condition <- condition
        all_diagnostics[[length(all_diagnostics) + 1]] <- diag_df
      }
    }
  }
}

combined_diagnostics <- rbindlist(all_diagnostics, fill = TRUE)
cat(sprintf("Total diagnostic records: %d\n\n", nrow(combined_diagnostics)))

# Extract relevant metrics and remove duplicates
diagnostics_summary <- combined_diagnostics %>%
  filter(y_type %in% c("original", "residual")) %>%
  select(embedding_name, n_sample, seed, condition, y_type, 
         r2_test, f_test_pvalue, globaltest_pvalue, dcor_pvalue,
         pcm_pvalue, rcot_pvalue, pcm_pvalue_split, rcot_pvalue_split) %>%
  as.data.table()

# Remove duplicates by taking the first row for each unique combination
diagnostics_summary <- diagnostics_summary[, .SD[1], 
                                          by = .(seed, n_sample, condition, embedding_name, y_type)]

# Reshape to have separate columns for original and residual R^2
r2_original <- diagnostics_summary[y_type == "original", 
                                   .(seed, n_sample, condition, embedding_name, r2_original = r2_test,
                                     pcm_pvalue, rcot_pvalue, pcm_pvalue_split, rcot_pvalue_split)]
r2_residual <- diagnostics_summary[y_type == "residual", 
                                   .(seed, n_sample, condition, embedding_name, r2_residual = r2_test,
                                     f_test_pvalue, globaltest_pvalue, dcor_pvalue)]

diagnostics_wide <- merge(r2_original, r2_residual, 
                          by = c("seed", "n_sample", "condition", "embedding_name"),
                          all = TRUE)

# ============================================================================
# Compute Per-Embedding Correlations
# ============================================================================
cat("=== Computing Per-Embedding Correlations ===\n")

correlation_results <- list()

for (cit_name in unique(combined_pvals$cit)) {
  for (emb in unique(combined_pvals$embedding)) {
    for (cond in conditions) {
      for (n_samp in n_samples) {
        # Get CIT p-values
        pval_subset <- combined_pvals[cit == cit_name & embedding == emb & 
                                      condition == cond & n_sample == n_samp]
        
        # Get diagnostics
        diag_subset <- diagnostics_wide[embedding_name == emb & condition == cond & 
                                       n_sample == n_samp]
        
        # Merge by seed
        merged_data <- merge(pval_subset, diag_subset, 
                            by = c("seed", "n_sample", "condition"),
                            all = FALSE)
        
        if (nrow(merged_data) < 3) next  # Need at least 3 observations
        
        # Compute correlations
        # Try globaltest correlation with error handling
        cor_globaltest <- tryCatch({
          cor(merged_data$p_value, merged_data$globaltest_pvalue, use = "complete.obs")
        }, error = function(e) {
          NA
        })
        
        cor_results <- data.frame(
          cit = cit_name,
          embedding = emb,
          condition = cond,
          n_sample = n_samp,
          cor_r2_orig = cor(merged_data$p_value, merged_data$r2_original, use = "complete.obs"),
          cor_r2_resid = cor(merged_data$p_value, merged_data$r2_residual, use = "complete.obs"),
          cor_ftest = cor(merged_data$p_value, merged_data$f_test_pvalue, use = "complete.obs"),
          cor_globaltest = cor_globaltest,
          cor_dcor = cor(merged_data$p_value, merged_data$dcor_pvalue, use = "complete.obs"),
          cor_pcm = cor(merged_data$p_value, merged_data$pcm_pvalue, use = "complete.obs"),
          cor_rcot = cor(merged_data$p_value, merged_data$rcot_pvalue, use = "complete.obs"),
          cor_pcm_split = cor(merged_data$p_value, merged_data$pcm_pvalue_split, use = "complete.obs"),
          cor_rcot_split = cor(merged_data$p_value, merged_data$rcot_pvalue_split, use = "complete.obs"),
          n_obs = nrow(merged_data)
        )
        
        correlation_results[[length(correlation_results) + 1]] <- cor_results
      }
    }
  }
}

correlation_df <- rbindlist(correlation_results)
cat(sprintf("Computed %d per-embedding correlations\n\n", nrow(correlation_df)))

# ============================================================================
# Compute Across-All-Embeddings Correlations
# ============================================================================
cat("=== Computing Across-All-Embeddings Correlations ===\n")

correlation_pooled_results <- list()

for (cit_name in unique(combined_pvals$cit)) {
  for (cond in conditions) {
    for (n_samp in n_samples) {
      # Get CIT p-values (all embeddings)
      pval_subset <- combined_pvals[cit == cit_name & condition == cond & n_sample == n_samp]
      
      # Get diagnostics (all embeddings)
      diag_subset <- diagnostics_wide[condition == cond & n_sample == n_samp]
      
      # Merge by seed and embedding
      merged_data <- merge(pval_subset, diag_subset, 
                          by.x = c("seed", "n_sample", "condition", "embedding"),
                          by.y = c("seed", "n_sample", "condition", "embedding_name"),
                          all = FALSE)
      
      if (nrow(merged_data) < 3) next
      
      # Compute correlations across all embeddings
      # Try globaltest correlation with error handling
      cor_globaltest <- tryCatch({
        cor(merged_data$p_value, merged_data$globaltest_pvalue, use = "complete.obs")
      }, error = function(e) {
        NA
      })
      
      cor_results <- data.frame(
        cit = cit_name,
        condition = cond,
        n_sample = n_samp,
        cor_r2_orig = cor(merged_data$p_value, merged_data$r2_original, use = "complete.obs"),
        cor_r2_resid = cor(merged_data$p_value, merged_data$r2_residual, use = "complete.obs"),
        cor_ftest = cor(merged_data$p_value, merged_data$f_test_pvalue, use = "complete.obs"),
        cor_globaltest = cor_globaltest,
        cor_dcor = cor(merged_data$p_value, merged_data$dcor_pvalue, use = "complete.obs"),
        cor_pcm = cor(merged_data$p_value, merged_data$pcm_pvalue, use = "complete.obs"),
        cor_rcot = cor(merged_data$p_value, merged_data$rcot_pvalue, use = "complete.obs"),
        cor_pcm_split = cor(merged_data$p_value, merged_data$pcm_pvalue_split, use = "complete.obs"),
        cor_rcot_split = cor(merged_data$p_value, merged_data$rcot_pvalue_split, use = "complete.obs"),
        n_obs = nrow(merged_data)
      )
      
      correlation_pooled_results[[length(correlation_pooled_results) + 1]] <- cor_results
    }
  }
}

correlation_pooled_df <- rbindlist(correlation_pooled_results)
cat(sprintf("Computed %d across-all-embeddings correlations\n\n", nrow(correlation_pooled_df)))

# ============================================================================
# Compute Seedwise Correlations Across Embeddings
# ============================================================================
cat("=== Computing Seedwise Correlations Across Embeddings ===\n")

correlation_seedwise_results <- list()

for (cit_name in unique(combined_pvals$cit)) {
  for (cond in conditions) {
    for (n_samp in n_samples) {
      for (s in seeds) {
        # Get CIT p-values for this seed (all embeddings)
        pval_subset <- combined_pvals[cit == cit_name & condition == cond & 
                                      n_sample == n_samp & seed == s]
        
        # Get diagnostics for this seed (all embeddings)
        diag_subset <- diagnostics_wide[condition == cond & n_sample == n_samp & seed == s]
        
        # Merge by embedding
        merged_data <- merge(pval_subset, diag_subset, 
                            by.x = c("seed", "n_sample", "condition", "embedding"),
                            by.y = c("seed", "n_sample", "condition", "embedding_name"),
                            all = FALSE)
        
        if (nrow(merged_data) < 3) next  # Need at least 3 embeddings
        
        # Compute correlations across embeddings for this seed
        cor_globaltest <- tryCatch({
          cor(merged_data$p_value, merged_data$globaltest_pvalue, use = "complete.obs")
        }, error = function(e) {
          NA
        })
        
        cor_results <- data.frame(
          cit = cit_name,
          condition = cond,
          n_sample = n_samp,
          seed = s,
          cor_r2_orig = cor(merged_data$p_value, merged_data$r2_original, use = "complete.obs"),
          cor_r2_resid = cor(merged_data$p_value, merged_data$r2_residual, use = "complete.obs"),
          cor_ftest = cor(merged_data$p_value, merged_data$f_test_pvalue, use = "complete.obs"),
          cor_globaltest = cor_globaltest,
          cor_dcor = cor(merged_data$p_value, merged_data$dcor_pvalue, use = "complete.obs"),
          cor_pcm = cor(merged_data$p_value, merged_data$pcm_pvalue, use = "complete.obs"),
          cor_rcot = cor(merged_data$p_value, merged_data$rcot_pvalue, use = "complete.obs"),
          cor_pcm_split = cor(merged_data$p_value, merged_data$pcm_pvalue_split, use = "complete.obs"),
          cor_rcot_split = cor(merged_data$p_value, merged_data$rcot_pvalue_split, use = "complete.obs"),
          n_obs = nrow(merged_data)
        )
        
        correlation_seedwise_results[[length(correlation_seedwise_results) + 1]] <- cor_results
      }
    }
  }
}

correlation_seedwise <- rbindlist(correlation_seedwise_results)
cat(sprintf("Computed %d seedwise correlations\n", nrow(correlation_seedwise)))

# Average across seeds
correlation_avg <- correlation_seedwise %>%
  group_by(cit, condition, n_sample) %>%
  summarize(
    avg_cor_r2_orig = mean(cor_r2_orig, na.rm = TRUE),
    avg_cor_r2_resid = mean(cor_r2_resid, na.rm = TRUE),
    avg_cor_ftest = mean(cor_ftest, na.rm = TRUE),
    avg_cor_globaltest = mean(cor_globaltest, na.rm = TRUE),
    avg_cor_dcor = mean(cor_dcor, na.rm = TRUE),
    avg_cor_pcm = mean(cor_pcm, na.rm = TRUE),
    avg_cor_rcot = mean(cor_rcot, na.rm = TRUE),
    avg_cor_pcm_split = mean(cor_pcm_split, na.rm = TRUE),
    avg_cor_rcot_split = mean(cor_rcot_split, na.rm = TRUE),
    n_seeds = sum(!is.na(cor_r2_orig)),
    .groups = "drop"
  )
cat(sprintf("Computed %d averaged correlations\n\n", nrow(correlation_avg)))

# ============================================================================
# Save Results
# ============================================================================
cat("=== Saving Results ===\n")

output_dir <- "/sc/home/marco.simnacher/dncitPaper/inst/diagnostic_embedding/results"

# Save per-embedding correlations
output_file_per_emb <- file.path(output_dir, "correlation_cit_diagnostics_per_embedding.csv")
fwrite(correlation_df, output_file_per_emb)
cat(sprintf("Saved per-embedding correlations: %s\n", output_file_per_emb))

# Save across-all-embeddings correlations
output_file_pooled <- file.path(output_dir, "correlation_cit_diagnostics_pooled.csv")
fwrite(correlation_pooled_df, output_file_pooled)
cat(sprintf("Saved pooled correlations: %s\n", output_file_pooled))

# Save avg-across-seeds correlations
output_file_avg <- file.path(output_dir, "correlation_cit_diagnostics_avg.csv")
fwrite(correlation_avg, output_file_avg)
cat(sprintf("Saved avg-across-seeds correlations: %s\n", output_file_avg))

cat("\n=== Summary ===\n")
cat("Per-embedding correlations:\n")
print(head(correlation_df))
cat("\nPooled correlations:\n")
print(head(correlation_pooled_df))

# ============================================================================
# Create Combined Diagnostic Plots (Boxplots + Correlation Lineplots)
# ============================================================================
cat("\n=== Creating Combined Diagnostic Plots ===\n")

# Define embeddings
constant_embeddings <- c("fastsurfer", "freesurfer", "condVAE", "medicalnet")
varying_embeddings <- c('medicalnet_ft', 'scratch')
available_embeddings <- unique(combined_pvals$embedding)

constant_embeddings_plot <- intersect(constant_embeddings, available_embeddings)
varying_embeddings_plot <- intersect(varying_embeddings, available_embeddings)

# Define embedding order: fastsurfer, freesurfer, condVAE, medicalnet, scratch, medicalnet_ft
embedding_order <- c("fastsurfer", "freesurfer", "condVAE", "medicalnet", "scratch", "medicalnet_ft")
all_embeddings_plot <- intersect(embedding_order, c(constant_embeddings_plot, varying_embeddings_plot))

# Create display name mapping for embeddings
embedding_display_names <- c(
  "fastsurfer" = "FAST",
  "freesurfer" = "Freesurfer",
  "condVAE" = "cVAE",
  "medicalnet" = "MedicalNet",
  "medicalnet_ft" = "MedicalNet-ft",
  "scratch" = "Scratch"
)

# Get display names for plotting
all_embeddings_display <- embedding_display_names[all_embeddings_plot]

# Create color palette using paletteer (same as pvalue_boxplots.R)
palet_discrete <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
color_palette <- setNames(
  palet_discrete[1:length(all_embeddings_plot)],
  all_embeddings_display
)

# Prepare diagnostic data for plotting (without CIT-specific merge)
# Use display names for embeddings and rename conditions
plot_data_diag <- diagnostics_wide %>%
  filter(embedding_name %in% all_embeddings_plot) %>%
  mutate(
    n_sample = factor(n_sample, levels = n_samples),
    condition = factor(condition, levels = conditions, labels = c("T1E", "Power")),
    embedding = factor(embedding_name, levels = all_embeddings_plot, labels = all_embeddings_display)
  )

# Debug: Check data availability
cat("\n=== Data Availability Check ===\n")
cat(sprintf("Total rows in plot_data_diag: %d\n", nrow(plot_data_diag)))
cat("\nRows by condition:\n")
print(table(plot_data_diag$condition))
cat("\nRows by condition and embedding:\n")
print(table(plot_data_diag$condition, plot_data_diag$embedding))
cat("\nNA counts for key metrics:\n")
cat(sprintf("  r2_original: %d NAs\n", sum(is.na(plot_data_diag$r2_original))))
cat(sprintf("  r2_residual: %d NAs\n", sum(is.na(plot_data_diag$r2_residual))))

# Prepare correlation data for lineplots
# Map CIT names for display and use display names for embeddings
correlation_df_plot <- correlation_df %>%
  filter(embedding %in% all_embeddings_plot) %>%
  mutate(
    n_sample = factor(n_sample, levels = n_samples),
    condition = factor(condition, levels = conditions, labels = c("T1E", "Power")),
    embedding = factor(embedding, levels = all_embeddings_plot, labels = all_embeddings_display),
    cit_label = case_when(
      grepl("RCOT", cit, ignore.case = TRUE) ~ "RCoT",
      grepl("comets_pcm", cit, ignore.case = TRUE) ~ "PCM",
      TRUE ~ as.character(cit)
    )
  )

# Create figures directory
figures_dir <- file.path(output_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Define diagnostic metrics to plot
# is_pvalue: TRUE for metrics that should use log10 scale
diagnostic_metrics <- list(
  r2_original = list(
    col = "r2_original",
    cor_col = "cor_r2_orig",
    y_label = expression(R^2 ~ "(Original)"),
    is_pvalue = FALSE
  ),
  r2_residual = list(
    col = "r2_residual",
    cor_col = "cor_r2_resid",
    y_label = expression(R^2 ~ "(Residual)"),
    is_pvalue = FALSE
  ),
  ftest = list(
    col = "f_test_pvalue",
    cor_col = "cor_ftest",
    y_label = "F-test P-value",
    is_pvalue = TRUE
  ),
  globaltest = list(
    col = "globaltest_pvalue",
    cor_col = "cor_globaltest",
    y_label = "Globaltest P-value",
    is_pvalue = TRUE
  ),
  dcor = list(
    col = "dcor_pvalue",
    cor_col = "cor_dcor",
    y_label = "Dcor P-value",
    is_pvalue = TRUE
  ),
  pcm = list(
    col = "pcm_pvalue",
    cor_col = "cor_pcm",
    y_label = "PCM P-value",
    is_pvalue = TRUE
  ),
  rcot = list(
    col = "rcot_pvalue",
    cor_col = "cor_rcot",
    y_label = "RCoT P-value",
    is_pvalue = TRUE
  ),
  pcm_split = list(
    col = "pcm_pvalue_split",
    cor_col = "cor_pcm_split",
    y_label = "PCM P-value (Split)",
    is_pvalue = TRUE
  ),
  rcot_split = list(
    col = "rcot_pvalue_split",
    cor_col = "cor_rcot_split",
    y_label = "RCoT P-value (Split)",
    is_pvalue = TRUE
  )
)

# Minimum p-value for log scale plotting
min_pval_plot <- 1e-16

# Common theme for all plots (matching pvalue_boxplots.R style)
base_theme <- theme_bw(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

# Create combined plots for each diagnostic metric
for (metric_name in names(diagnostic_metrics)) {
  cat(sprintf("\nCreating combined plot for: %s\n", metric_name))
  
  metric_info <- diagnostic_metrics[[metric_name]]
  
  tryCatch({
    # ===== TOP ROW: Boxplots =====
    boxplot_data <- plot_data_diag %>%
      select(seed, n_sample, condition, embedding, !!sym(metric_info$col)) %>%
      rename(metric_value = !!sym(metric_info$col)) %>%
      filter(!is.na(metric_value))
    
    # For p-value metrics, apply floor at min_pval_plot
    if (metric_info$is_pvalue) {
      boxplot_data <- boxplot_data %>%
        mutate(metric_value = ifelse(metric_value < min_pval_plot | metric_value == 0, 
                                      min_pval_plot, metric_value))
    }
    
    if (nrow(boxplot_data) == 0) {
      cat(sprintf("  Skipping %s - no data available\n", metric_name))
      next
    }
    
    p_boxplot <- ggplot(boxplot_data, aes(x = n_sample, y = metric_value, fill = embedding)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
      scale_y_log10()
      facet_wrap(~ condition, ncol = 2) +
      labs(
        y = metric_info$y_label,
        fill = "Embedding"
      ) +
      base_theme +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
      ) +
      scale_fill_manual(values = color_palette)
    
    # Add log10 scale for p-value metrics
    if (metric_info$is_pvalue) {
      p_boxplot <- p_boxplot +
        scale_y_log10(
          breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 1),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          limits = c(min_pval_plot, 1)
        ) +
        geom_hline(yintercept = 0.05, color = "black", alpha = 0.6)
    }
    
    # ===== MIDDLE ROW: PCM Correlation Lineplots =====
    # Get correlation data for this metric - PCM only
    cor_data_pcm <- correlation_df_plot %>%
      filter(cit_label == "PCM") %>%
      select(cit, cit_label, embedding, condition, n_sample, !!sym(metric_info$cor_col)) %>%
      rename(correlation = !!sym(metric_info$cor_col)) %>%
      filter(!is.na(correlation))
    
    # Extract averaged correlation for PCM
    avg_cor_col_name <- paste0("avg_", metric_info$cor_col)
    avg_cor_pcm <- correlation_avg %>%
      filter(grepl("comets_pcm", cit, ignore.case = TRUE)) %>%
      mutate(
        n_sample = factor(n_sample, levels = n_samples),
        condition = factor(condition, levels = conditions, labels = c("T1E", "Power"))
      ) %>%
      rename(avg_correlation = !!sym(avg_cor_col_name))
    
    # ===== BOTTOM ROW: RCoT Correlation Lineplots =====
    # Get correlation data for this metric - RCoT only
    cor_data_rcot <- correlation_df_plot %>%
      filter(cit_label == "RCoT") %>%
      select(cit, cit_label, embedding, condition, n_sample, !!sym(metric_info$cor_col)) %>%
      rename(correlation = !!sym(metric_info$cor_col)) %>%
      filter(!is.na(correlation))
    
    # Extract averaged correlation for RCoT
    avg_cor_rcot <- correlation_avg %>%
      filter(grepl("RCOT", cit, ignore.case = TRUE)) %>%
      mutate(
        n_sample = factor(n_sample, levels = n_samples),
        condition = factor(condition, levels = conditions, labels = c("T1E", "Power"))
      ) %>%
      rename(avg_correlation = !!sym(avg_cor_col_name))
    
    if (nrow(cor_data_pcm) == 0 && nrow(cor_data_rcot) == 0) {
      cat(sprintf("  No correlation data for %s - creating boxplot only\n", metric_name))
      combined_plot <- p_boxplot
    } else {
      # Combine PCM and RCoT correlation data with CIT label for unified legend
      cor_data_combined <- bind_rows(
        cor_data_pcm %>% mutate(CIT = "PCM"),
        cor_data_rcot %>% mutate(CIT = "RCoT")
      ) %>%
        mutate(CIT = factor(CIT, levels = c("PCM", "RCoT")))
      
      # PCM correlation lineplot (middle row)
      p_lineplot_pcm <- ggplot(cor_data_pcm %>% mutate(CIT = "PCM"), 
                                aes(x = n_sample, y = correlation, 
                                    color = embedding, group = embedding,
                                    linetype = CIT, shape = CIT)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2.5) +
        geom_line(data = avg_cor_pcm, aes(x = n_sample, y = avg_correlation, group = 1),
                  inherit.aes = FALSE, color = "black", linewidth = 1.2) +
        geom_point(data = avg_cor_pcm, aes(x = n_sample, y = avg_correlation),
                   inherit.aes = FALSE, color = "black", size = 3, shape = 18) +
        facet_wrap(~ condition, ncol = 2) +
        labs(
          y = "Correlation (PCM)",
          color = "Embedding",
          linetype = "CIT",
          shape = "CIT"
        ) +
        base_theme +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right"
        ) +
        scale_color_manual(values = color_palette) +
        scale_linetype_manual(values = c("PCM" = "solid", "RCoT" = "dashed")) +
        scale_shape_manual(values = c("PCM" = 16, "RCoT" = 17)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.7)
      
      # RCoT correlation lineplot (bottom row)
      p_lineplot_rcot <- ggplot(cor_data_rcot %>% mutate(CIT = "RCoT"), 
                                 aes(x = n_sample, y = correlation, 
                                     color = embedding, group = embedding,
                                     linetype = CIT, shape = CIT)) +
        geom_line(linewidth = 0.8) +
        geom_point(size = 2.5) +
        geom_line(data = avg_cor_rcot, aes(x = n_sample, y = avg_correlation, group = 1),
                  inherit.aes = FALSE, color = "black", linewidth = 1.2) +
        geom_point(data = avg_cor_rcot, aes(x = n_sample, y = avg_correlation),
                   inherit.aes = FALSE, color = "black", size = 3, shape = 18) +
        facet_wrap(~ condition, ncol = 2) +
        labs(
          x = "Sample Size",
          y = "Correlation (RCoT)",
          color = "Embedding",
          linetype = "CIT",
          shape = "CIT"
        ) +
        base_theme +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right"
        ) +
        scale_color_manual(values = color_palette) +
        scale_linetype_manual(values = c("PCM" = "solid", "RCoT" = "dashed")) +
        scale_shape_manual(values = c("PCM" = 16, "RCoT" = 17)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.7)
      
      # ===== Combine plots using patchwork (3 rows) =====
      # Use guides = "collect" to merge legends, then override legend aesthetics
      combined_plot <- p_boxplot / p_lineplot_pcm / p_lineplot_rcot +
        plot_layout(heights = c(1, 1, 1), guides = "collect") &
        theme(legend.position = "right") &
        guides(
          fill = guide_legend(order = 1, title = "Embedding"),
          color = guide_legend(order = 1, title = "Embedding"),
          linetype = guide_legend(order = 2, title = "CIT"),
          shape = guide_legend(order = 2, title = "CIT")
        )
    }
    
    # Save combined plot
    png_path <- file.path(figures_dir, sprintf("diagnostic_%s_%d_%d.png", 
                                               metric_name, min(seeds), max(seeds)))
    ggsave(png_path, plot = combined_plot, width = 14, height = 14, units = "in", dpi = 300)
    cat(sprintf("  Saved: %s\n", basename(png_path)))
    
    pdf_path <- file.path(figures_dir, sprintf("diagnostic_%s_%d_%d.pdf", 
                                               metric_name, min(seeds), max(seeds)))
    ggsave(pdf_path, plot = combined_plot, width = 14, height = 14, units = "in")
    
  }, error = function(e) {
    cat(sprintf("  Warning: Plot failed for %s: %s\n", metric_name, e$message))
  })
}

cat("\n=== Diagnostic Visualization Complete ===\n")
cat(sprintf("All plots saved to: %s\n", figures_dir))

# ============================================================================
# Create Rank Scatter Plots (New Separate Figures)
# ============================================================================
cat("\n=== Creating Rank Scatter Plots ===\n")

# Merge diagnostic data with CIT p-values for rank computation
# Need to merge diagnostics_wide with combined_pvals by seed, n_sample, condition, embedding

# Prepare CIT p-values with display names for embeddings
cit_pvals_for_ranks <- combined_pvals %>%
  filter(embedding %in% all_embeddings_plot) %>%
  mutate(
    n_sample = factor(n_sample, levels = n_samples),
    condition = factor(condition, levels = conditions, labels = c("T1E", "Power")),
    embedding = factor(embedding, levels = all_embeddings_plot, labels = all_embeddings_display),
    cit_label = case_when(
      grepl("RCOT", cit, ignore.case = TRUE) ~ "RCoT",
      grepl("comets_pcm", cit, ignore.case = TRUE) ~ "PCM",
      TRUE ~ as.character(cit)
    )
  )

# Merge with diagnostic data
rank_data_base <- plot_data_diag %>%
  inner_join(
    cit_pvals_for_ranks %>% select(seed, n_sample, condition, embedding, cit_label, p_value),
    by = c("seed", "n_sample", "condition", "embedding")
  )

# Create rank scatter plots for each diagnostic metric
for (metric_name in names(diagnostic_metrics)) {
  cat(sprintf("\nCreating rank scatter plot for: %s\n", metric_name))
  
  metric_info <- diagnostic_metrics[[metric_name]]
  
  tryCatch({
    # Prepare data for this metric
    metric_col <- metric_info$col
    
    rank_data <- rank_data_base %>%
      select(seed, n_sample, condition, embedding, cit_label, 
             metric_value = !!sym(metric_col), p_value) %>%
      filter(!is.na(metric_value), !is.na(p_value))
    
    if (nrow(rank_data) == 0) {
      cat(sprintf("  Skipping %s - no data available for ranks\n", metric_name))
      next
    }
    
    # Compute ranks within each (seed, n_sample, condition, cit_label) group
    # Rank across embeddings (1 to 6)
    rank_data <- rank_data %>%
      group_by(seed, n_sample, condition, cit_label) %>%
      mutate(
        metric_rank = rank(metric_value, ties.method = "average"),
        pvalue_rank = rank(p_value, ties.method = "average")
      ) %>%
      ungroup()
    
    # Calculate Spearman correlation for each (condition, cit_label) combination
    spearman_results <- rank_data %>%
      group_by(condition, cit_label) %>%
      summarize(
        spearman_rho = cor(metric_rank, pvalue_rank, method = "spearman", use = "complete.obs"),
        spearman_test = list(cor.test(metric_rank, pvalue_rank, method = "spearman")),
        .groups = "drop"
      ) %>%
      mutate(
        spearman_pval = sapply(spearman_test, function(x) x$p.value),
        annotation = sprintf("rho = %.3f\np = %.2e", spearman_rho, spearman_pval)
      )
    
    # ===== TOP ROW: Boxplots (same as before) =====
    boxplot_data_ranks <- plot_data_diag %>%
      select(seed, n_sample, condition, embedding, !!sym(metric_col)) %>%
      rename(metric_value = !!sym(metric_col)) %>%
      filter(!is.na(metric_value))
    
    if (metric_info$is_pvalue) {
      boxplot_data_ranks <- boxplot_data_ranks %>%
        mutate(metric_value = ifelse(metric_value < min_pval_plot | metric_value == 0, 
                                      min_pval_plot, metric_value))
    }
    
    p_boxplot_ranks <- ggplot(boxplot_data_ranks, aes(x = n_sample, y = metric_value, fill = embedding)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
      facet_wrap(~ condition, ncol = 2) +
      labs(
        y = metric_info$y_label,
        fill = "Embedding"
      ) +
      base_theme +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
      ) +
      scale_fill_manual(values = color_palette)
    
    if (metric_info$is_pvalue) {
      p_boxplot_ranks <- p_boxplot_ranks +
        scale_y_log10(
          breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 1),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          limits = c(min_pval_plot, 1)
        ) +
        geom_hline(yintercept = 0.05, color = "black", alpha = 0.6)
    }
    
    # ===== MIDDLE ROW: PCM Rank Scatter =====
    rank_data_pcm <- rank_data %>% filter(cit_label == "PCM")
    spearman_pcm <- spearman_results %>% filter(cit_label == "PCM")
    
    p_scatter_pcm <- ggplot(rank_data_pcm, aes(x = metric_rank, y = pvalue_rank, color = embedding)) +
      geom_jitter(alpha = 0.5, width = 0.1, height = 0.1, size = 1.5) +
      geom_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black", linetype = "dashed", linewidth = 0.8) +
      geom_text(data = spearman_pcm, aes(x = Inf, y = Inf, label = annotation),
                inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
      facet_wrap(~ condition, ncol = 2) +
      labs(
        y = "Rank of PCM P-value",
        color = "Embedding"
      ) +
      base_theme +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
      ) +
      scale_color_manual(values = color_palette) +
      scale_x_continuous(breaks = 1:6, limits = c(0.5, 6.5)) +
      scale_y_continuous(breaks = 1:6, limits = c(0.5, 6.5))
    
    # ===== BOTTOM ROW: RCoT Rank Scatter =====
    rank_data_rcot <- rank_data %>% filter(cit_label == "RCoT")
    spearman_rcot <- spearman_results %>% filter(cit_label == "RCoT")
    
    p_scatter_rcot <- ggplot(rank_data_rcot, aes(x = metric_rank, y = pvalue_rank, color = embedding)) +
      geom_jitter(alpha = 0.5, width = 0.1, height = 0.1, size = 1.5) +
      geom_smooth(method = "lm", se = FALSE, aes(group = 1), color = "black", linetype = "dashed", linewidth = 0.8) +
      geom_text(data = spearman_rcot, aes(x = Inf, y = Inf, label = annotation),
                inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
      facet_wrap(~ condition, ncol = 2) +
      labs(
        x = "Rank of Diagnostic Metric",
        y = "Rank of RCoT P-value",
        color = "Embedding"
      ) +
      base_theme +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "right"
      ) +
      scale_color_manual(values = color_palette) +
      scale_x_continuous(breaks = 1:6, limits = c(0.5, 6.5)) +
      scale_y_continuous(breaks = 1:6, limits = c(0.5, 6.5))
    
    # ===== Combine rank scatter plots using patchwork (3 rows) =====
    combined_rank_plot <- p_boxplot_ranks / p_scatter_pcm / p_scatter_rcot +
      plot_layout(heights = c(1, 1, 1), guides = "collect") &
      theme(legend.position = "right")
    
    # Save rank scatter plot
    png_path_ranks <- file.path(figures_dir, sprintf("diagnostic_%s_ranks_%d_%d.png", 
                                                      metric_name, min(seeds), max(seeds)))
    ggsave(png_path_ranks, plot = combined_rank_plot, width = 14, height = 14, units = "in", dpi = 300)
    cat(sprintf("  Saved: %s\n", basename(png_path_ranks)))
    
    pdf_path_ranks <- file.path(figures_dir, sprintf("diagnostic_%s_ranks_%d_%d.pdf", 
                                                      metric_name, min(seeds), max(seeds)))
    ggsave(pdf_path_ranks, plot = combined_rank_plot, width = 14, height = 14, units = "in")
    
  }, error = function(e) {
    cat(sprintf("  Warning: Rank scatter plot failed for %s: %s\n", metric_name, e$message))
  })
}

cat("\n=== Rank Scatter Plots Complete ===\n")

# ============================================================================
# Create Rank Correlation Lineplot Figures (Spearman correlation per sample size)
# ============================================================================
cat("\n=== Creating Rank Correlation Lineplot Figures ===\n")

# Create rank correlation lineplots for each diagnostic metric
for (metric_name in names(diagnostic_metrics)) {
  cat(sprintf("\nCreating rank correlation lineplot for: %s\n", metric_name))
  
  metric_info <- diagnostic_metrics[[metric_name]]
  
  tryCatch({
    # Prepare data for this metric
    metric_col <- metric_info$col
    
    rank_data <- rank_data_base %>%
      select(seed, n_sample, condition, embedding, cit_label, 
             metric_value = !!sym(metric_col), p_value) %>%
      filter(!is.na(metric_value), !is.na(p_value))
    
    if (nrow(rank_data) == 0) {
      cat(sprintf("  Skipping %s - no data available for rank correlations\n", metric_name))
      next
    }
    
    # Compute ranks within each (seed, n_sample, condition, cit_label) group
    rank_data <- rank_data %>%
      group_by(seed, n_sample, condition, cit_label) %>%
      mutate(
        metric_rank = rank(metric_value, ties.method = "average"),
        pvalue_rank = rank(p_value, ties.method = "average")
      ) %>%
      ungroup()
    
    # Compute Spearman correlation per (condition, n_sample, cit_label, embedding)
    # Pool all seeds together for each combination
    spearman_per_sample <- rank_data %>%
      group_by(condition, n_sample, cit_label, embedding) %>%
      summarize(
        spearman_rho = cor(metric_rank, pvalue_rank, method = "spearman", use = "complete.obs"),
        n_obs = n(),
        .groups = "drop"
      )
    
    # Compute correlation seedwise across embeddings, then summarize across seeds
    # (This reflects "across-embedding" rank association per sample size.)
    spearman_seedwise <- rank_data %>%
      group_by(seed, condition, n_sample, cit_label) %>%
      summarize(
        spearman_rho = suppressWarnings(cor(metric_rank, pvalue_rank, method = "spearman", use = "complete.obs")),
        .groups = "drop"
      )

    spearman_avg <- spearman_seedwise %>%
      group_by(condition, n_sample, cit_label) %>%
      summarize(
        avg_spearman_rho = mean(spearman_rho, na.rm = TRUE),
        n_seeds = sum(!is.na(spearman_rho)),
        .groups = "drop"
      )
    
    # ===== TOP ROW: Boxplots (same as before) =====
    boxplot_data_rankcor <- plot_data_diag %>%
      select(seed, n_sample, condition, embedding, !!sym(metric_col)) %>%
      rename(metric_value = !!sym(metric_col)) %>%
      filter(!is.na(metric_value))
    
    if (metric_info$is_pvalue) {
      boxplot_data_rankcor <- boxplot_data_rankcor %>%
        mutate(metric_value = ifelse(metric_value < min_pval_plot | metric_value == 0, 
                                      min_pval_plot, metric_value))
    }
    
    p_boxplot_rankcor <- ggplot(boxplot_data_rankcor, aes(x = n_sample, y = metric_value, fill = embedding)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
      facet_wrap(~ condition, ncol = 2) +
      labs(
        y = metric_info$y_label,
        fill = "Embedding"
      ) +
      base_theme +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
      ) +
      scale_fill_manual(values = color_palette)
    
    if (metric_info$is_pvalue) {
      p_boxplot_rankcor <- p_boxplot_rankcor +
        scale_y_log10(
          breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 1),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          limits = c(min_pval_plot, 1)
        ) +
        geom_hline(yintercept = 0.05, color = "black", alpha = 0.6)
    }
    
    # ===== MIDDLE ROW: PCM Spearman Correlation Lineplot =====
    spearman_pcm <- spearman_per_sample %>% 
      filter(cit_label == "PCM") %>%
      mutate(CIT = "PCM")
    
    spearman_avg_pcm <- spearman_avg %>% filter(cit_label == "PCM")
    
    p_rankcor_pcm <- ggplot(spearman_pcm, aes(x = n_sample, y = spearman_rho, 
                                               color = embedding, group = embedding,
                                               linetype = CIT, shape = CIT)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2.5) +
      geom_line(data = spearman_avg_pcm, aes(x = n_sample, y = avg_spearman_rho, group = 1),
                inherit.aes = FALSE, color = "black", linewidth = 1.2) +
      geom_point(data = spearman_avg_pcm, aes(x = n_sample, y = avg_spearman_rho),
                 inherit.aes = FALSE, color = "black", size = 3, shape = 18) +
      facet_wrap(~ condition, ncol = 2) +
      labs(
        y = "Spearman Correlation (PCM)",
        color = "Embedding",
        linetype = "CIT",
        shape = "CIT"
      ) +
      base_theme +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right"
      ) +
      scale_color_manual(values = color_palette) +
      scale_linetype_manual(values = c("PCM" = "solid", "RCoT" = "dashed")) +
      scale_shape_manual(values = c("PCM" = 16, "RCoT" = 17)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.7) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25))
    
    # ===== BOTTOM ROW: RCoT Spearman Correlation Lineplot =====
    spearman_rcot <- spearman_per_sample %>% 
      filter(cit_label == "RCoT") %>%
      mutate(CIT = "RCoT")
    
    spearman_avg_rcot <- spearman_avg %>% filter(cit_label == "RCoT")
    
    p_rankcor_rcot <- ggplot(spearman_rcot, aes(x = n_sample, y = spearman_rho, 
                                                 color = embedding, group = embedding,
                                                 linetype = CIT, shape = CIT)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2.5) +
      geom_line(data = spearman_avg_rcot, aes(x = n_sample, y = avg_spearman_rho, group = 1),
                inherit.aes = FALSE, color = "black", linewidth = 1.2) +
      geom_point(data = spearman_avg_rcot, aes(x = n_sample, y = avg_spearman_rho),
                 inherit.aes = FALSE, color = "black", size = 3, shape = 18) +
      facet_wrap(~ condition, ncol = 2) +
      labs(
        x = "Sample Size",
        y = "Spearman Correlation (RCoT)",
        color = "Embedding",
        linetype = "CIT",
        shape = "CIT"
      ) +
      base_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      scale_color_manual(values = color_palette) +
      scale_linetype_manual(values = c("PCM" = "solid", "RCoT" = "dashed")) +
      scale_shape_manual(values = c("PCM" = 16, "RCoT" = 17)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", alpha = 0.7) +
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25))
    
    # ===== Combine rank correlation lineplot using patchwork (3 rows) =====
    combined_rankcor_plot <- p_boxplot_rankcor / p_rankcor_pcm / p_rankcor_rcot +
      plot_layout(heights = c(1, 1, 1), guides = "collect") &
      theme(legend.position = "right") &
      guides(
        fill = guide_legend(order = 1, title = "Embedding"),
        color = guide_legend(order = 1, title = "Embedding"),
        linetype = guide_legend(order = 2, title = "CIT"),
        shape = guide_legend(order = 2, title = "CIT")
      )
    
    # Save rank correlation lineplot
    png_path_rankcor <- file.path(figures_dir, sprintf("diagnostic_%s_rankcor_%d_%d.png", 
                                                        metric_name, min(seeds), max(seeds)))
    ggsave(png_path_rankcor, plot = combined_rankcor_plot, width = 14, height = 14, units = "in", dpi = 300)
    cat(sprintf("  Saved: %s\n", basename(png_path_rankcor)))
    
    pdf_path_rankcor <- file.path(figures_dir, sprintf("diagnostic_%s_rankcor_%d_%d.pdf", 
                                                        metric_name, min(seeds), max(seeds)))
    ggsave(pdf_path_rankcor, plot = combined_rankcor_plot, width = 14, height = 14, units = "in")
    
  }, error = function(e) {
    cat(sprintf("  Warning: Rank correlation lineplot failed for %s: %s\n", metric_name, e$message))
  })
}

cat("\n=== Rank Correlation Lineplots Complete ===\n")

# ============================================================================
# NEW FIGURE 1: ECDF diagnostics (colored by embedding; linetype by n for scratch/ft)
# ============================================================================
cat("\n=== Creating ECDF Diagnostic Plots ===\n")

library(scales)

# which embeddings vary with sample size (training-size dependence)
varying_embeddings_disp <- c("Scratch", "MedicalNet-ft")
fixed_embeddings_disp   <- setdiff(levels(plot_data_diag$embedding), varying_embeddings_disp)

# For fixed embeddings: use only one n_sample to avoid repeated identical ECDFs
n_ref_for_fixed <- max(n_samples)  # you can change this to e.g. min(n_samples)

# Linetype mapping for varying embeddings: different per n_sample
linetype_map <- setNames(
  c("dashed", "dotted", "solid")[seq_along(n_samples)],
  as.character(n_samples)
)

# Function to create and save ECDF plot for a given diagnostic metric
make_ecdf_diag_plot <- function(metric_name, metric_info) {
  metric_col <- metric_info$col

  ecdf_data <- plot_data_diag %>%
    select(seed, n_sample, condition, embedding, metric_value = !!sym(metric_col)) %>%
    filter(!is.na(metric_value)) %>%
    # For fixed embeddings keep only one reference n_sample
    filter(
      (as.character(embedding) %in% varying_embeddings_disp) |
        (as.character(embedding) %in% fixed_embeddings_disp & as.character(n_sample) == as.character(n_ref_for_fixed))
    ) %>%
    mutate(
      # Only scratch/mednet-ft have linetype by n_sample; others are fixed (solid)
      linetype_by_n = ifelse(as.character(embedding) %in% varying_embeddings_disp,
                             as.character(n_sample),
                             "fixed")
    )

  if (nrow(ecdf_data) == 0) {
    message(sprintf("ECDF plot skipped for %s (no data)", metric_name))
    return(NULL)
  }

  # Floor p-values for log plotting if needed
  if (isTRUE(metric_info$is_pvalue)) {
    ecdf_data <- ecdf_data %>%
      mutate(metric_value = ifelse(metric_value < min_pval_plot | metric_value == 0,
                                   min_pval_plot, metric_value))
  }

  # linetype scale: fixed = solid, varying = mapped by n_sample
  lt_vals <- c("fixed" = "solid", linetype_map)

  p <- ggplot(ecdf_data, aes(x = metric_value, color = embedding, linetype = linetype_by_n)) +
    stat_ecdf(geom = "step", linewidth = 0.9, alpha = 0.95) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = metric_info$y_label,
      y = "ECDF",
      color = "Embedding",
      linetype = "Train n (Scratch / MedNet-ft)"
    ) +
    base_theme +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16)
    ) +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = lt_vals)

  # If p-values: log scale + reference line at 0.05 (optional)
  if (FALSE){#isTRUE(metric_info$is_pvalue)) {
    p <- p +
      scale_x_log10(
        breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 1),
        limits = c(min_pval_plot, 1)
      ) +
      geom_vline(xintercept = 0.05, color = "black", alpha = 0.4)
  }

  # Save
  pdf_path <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_%d_%d_eps%.1f.pdf", metric_name, min(seeds), max(seeds), eps_sigmaY))
  png_path <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_%d_%d_eps%.1f.png", metric_name, min(seeds), max(seeds), eps_sigmaY))
  ggsave(pdf_path, p, width = 14, height = 6, units = "in")
  ggsave(png_path, p, width = 14, height = 6, units = "in", dpi = 300)

  message(sprintf("Saved ECDF diagnostic plot: %s", basename(pdf_path)))
  p
}

# Create ECDF plots for all diagnostic metrics (or subset here)
for (metric_name in names(diagnostic_metrics)) {
  metric_info <- diagnostic_metrics[[metric_name]]
  tryCatch(
    make_ecdf_diag_plot(metric_name, metric_info),
    error = function(e) message(sprintf("ECDF failed for %s: %s", metric_name, e$message))
  )
}

cat("\n=== ECDF Diagnostic Plots Complete ===\n")

# ============================================================================
# NEW FIGURE 1 variants:
#   A) power_log10: log10 x-scale ONLY for Power panel (p-values)
#   B) zoom_0_0.1: keep full plot + add Power zoom panel [0, 0.1]
# ============================================================================
cat("\n=== Creating ECDF Diagnostic Plots (with variants) ===\n")

library(scales)
library(patchwork)

# which embeddings vary with sample size (training-size dependence)
varying_embeddings_disp <- c("Scratch", "MedicalNet-ft")
fixed_embeddings_disp   <- setdiff(levels(plot_data_diag$embedding), varying_embeddings_disp)

# For fixed embeddings: use only one n_sample to avoid repeated identical ECDFs
n_ref_for_fixed <- max(n_samples)

# Linetype mapping for varying embeddings: different per n_sample
linetype_map <- setNames(
  c("dashed", "dotted", "solid")[seq_along(n_samples)],
  as.character(n_samples)
)

# helper: build base ECDF layer for a given data subset
.build_ecdf_plot <- function(dat, metric_info, lt_vals, add_legend = TRUE) {
  p <- ggplot(dat, aes(x = metric_value, color = embedding, linetype = linetype_by_n)) +
    stat_ecdf(geom = "step", linewidth = 0.9, alpha = 0.95) +
    labs(
      x = metric_info$y_label,
      y = "ECDF",
      color = "Embedding",
      linetype = "Train n (Scratch / MedNet-ft)"
    ) +
    base_theme +
    theme(
      legend.position = if (add_legend) "right" else "none",
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16)
    ) +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = lt_vals)

  # Optional alpha threshold reference (works for both linear/log scales)
  if (isTRUE(metric_info$is_pvalue)) {
    p <- p + geom_vline(xintercept = 0.05, color = "black", alpha = 0.4, linetype = "dashed")
  }
  p
}

# Function to create and save ECDF plots for a given diagnostic metric
make_ecdf_diag_plot <- function(metric_name, metric_info,
                                make_power_log10 = TRUE,
                                make_zoom_0_0.1 = TRUE) {

  metric_col <- metric_info$col

  ecdf_data <- plot_data_diag %>%
    select(seed, n_sample, condition, embedding, metric_value = !!sym(metric_col)) %>%
    filter(!is.na(metric_value)) %>%
    # For fixed embeddings keep only one reference n_sample
    filter(
      (as.character(embedding) %in% varying_embeddings_disp) |
        (as.character(embedding) %in% fixed_embeddings_disp &
           as.character(n_sample) == as.character(n_ref_for_fixed))
    ) %>%
    mutate(
      linetype_by_n = ifelse(as.character(embedding) %in% varying_embeddings_disp,
                             as.character(n_sample),
                             "fixed")
    )

  if (nrow(ecdf_data) == 0) {
    message(sprintf("ECDF plot skipped for %s (no data)", metric_name))
    return(NULL)
  }

  # Floor p-values (important for log scale)
  if (isTRUE(metric_info$is_pvalue)) {
    ecdf_data <- ecdf_data %>%
      mutate(metric_value = ifelse(metric_value < min_pval_plot | metric_value == 0,
                                   min_pval_plot, metric_value))
  }

  # linetype scale: fixed = solid, varying = mapped by n_sample
  lt_vals <- c("fixed" = "solid", linetype_map)

  # -----------------------
  # BASE ECDF (your current one): facet_wrap, linear x
  # -----------------------
  p_base <- ggplot(ecdf_data, aes(x = metric_value, color = embedding, linetype = linetype_by_n)) +
    stat_ecdf(geom = "step", linewidth = 0.9, alpha = 0.95) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = metric_info$y_label,
      y = "ECDF",
      color = "Embedding",
      linetype = "Train n (Scratch / MedNet-ft)"
    ) +
    base_theme +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16)
    ) +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = lt_vals)

  if (isTRUE(metric_info$is_pvalue)) {
    p_base <- p_base + geom_vline(xintercept = 0.05, color = "black", alpha = 0.4, linetype = "dashed")
  }

  # Save base
  pdf_path <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_%d_%d_eps%.1f.pdf", metric_name, min(seeds), max(seeds), eps_sigmaY))
  png_path <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_%d_%d_eps%.1f.png", metric_name, min(seeds), max(seeds), eps_sigmaY))
  ggsave(pdf_path, p_base, width = 14, height = 6, units = "in")
  ggsave(png_path, p_base, width = 14, height = 6, units = "in", dpi = 300)
  message(sprintf("Saved ECDF (base): %s", basename(pdf_path)))

  # -----------------------
  # VARIANT (i): Power log10 x-scale only
  # Implemented as two separate panels combined with patchwork:
  #   left: T1E linear [0,1]
  #   right: Power log10 (p-values only); linear otherwise
  # -----------------------
  if (make_power_log10) {
    dat_t1e   <- ecdf_data %>% filter(condition == "T1E")
    dat_power <- ecdf_data %>% filter(condition == "Power")

    p_t1e <- .build_ecdf_plot(dat_t1e, metric_info, lt_vals, add_legend = FALSE) +
      labs(title = NULL) +
      theme(strip.background = element_blank(),
            strip.text = element_blank()) +
      ggtitle("T1E")

    p_power <- .build_ecdf_plot(dat_power, metric_info, lt_vals, add_legend = TRUE) +
      labs(title = NULL) +
      theme(strip.background = element_blank(),
            strip.text = element_blank()) +
      ggtitle("Power")

    # apply scaling ONLY to power panel if p-value metric
    if (isTRUE(metric_info$is_pvalue)) {
      p_power <- p_power +
        scale_x_log10(
          breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 1),
          limits = c(min_pval_plot, 1)
        )
      # keep T1E linear (default)
    }

    p_powerlog <- (p_t1e | p_power) +
      plot_layout(guides = "collect") &
      theme(legend.position = "right")

    pdf_path_pl <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_powerlog_%d_%d_eps%.1f.pdf", metric_name, min(seeds), max(seeds), eps_sigmaY))
    png_path_pl <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_powerlog_%d_%d_eps%.1f.png", metric_name, min(seeds), max(seeds), eps_sigmaY))
    ggsave(pdf_path_pl, p_powerlog, width = 14, height = 6, units = "in")
    ggsave(png_path_pl, p_powerlog, width = 14, height = 6, units = "in", dpi = 300)
    message(sprintf("Saved ECDF (power log10): %s", basename(pdf_path_pl)))
  }

  # -----------------------
  # VARIANT (ii): Zoom into [0, 0.1] while keeping full plot
  # Implemented as: full (both facets) + extra Power zoom panel on the right
  # -----------------------
  if (make_zoom_0_0.1) {
    dat_power <- ecdf_data %>% filter(condition == "Power")

    # Power-only zoom panel: same ECDFs, but xlim [0,0.1]
    p_power_zoom <- .build_ecdf_plot(dat_power, metric_info, lt_vals, add_legend = FALSE) +
      #coord_cartesian(xlim = c(0, 0.1)) +
      scale_x_log10(limits= c(min_pval_plot, 0.1)) +
      ggtitle("Power (zoom: x in [0, 0.1])") +
      theme(plot.title = element_text(size = 16))

    # Combine: keep original full two-facet plot + zoom
    p_zoom_combo <- (p_base | p_power_zoom) +
      plot_layout(widths = c(2.2, 1), guides = "collect") &
      theme(legend.position = "right")

    pdf_path_z <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_zoom01_%d_%d_eps%.1f.pdf", metric_name, min(seeds), max(seeds), eps_sigmaY))
    png_path_z <- file.path(figures_dir, sprintf("diagnostic_ecdf_%s_zoom01_%d_%d_eps%.1f.png", metric_name, min(seeds), max(seeds), eps_sigmaY))
    ggsave(pdf_path_z, p_zoom_combo, width = 18, height = 6, units = "in")
    ggsave(png_path_z, p_zoom_combo, width = 18, height = 6, units = "in", dpi = 300)
    message(sprintf("Saved ECDF (zoom [0,0.1] + full): %s", basename(pdf_path_z)))
  }

  invisible(TRUE)
}

# Create ECDF plots (base + variants) for all diagnostic metrics
for (metric_name in names(diagnostic_metrics)) {
  metric_info <- diagnostic_metrics[[metric_name]]
  tryCatch(
    make_ecdf_diag_plot(metric_name, metric_info,
                        make_power_log10 = TRUE,
                        make_zoom_0_0.1 = TRUE),
    error = function(e) message(sprintf("ECDF failed for %s: %s", metric_name, e$message))
  )
}

cat("\n=== ECDF Diagnostic Plots (with variants) Complete ===\n")

library(patchwork)

# ----------------------------
# helper: get base metric name from split metric name
# ----------------------------
.get_base_metric <- function(metric_name) {
  if (grepl("_split$", metric_name)) {
    sub("_split$", "", metric_name)
  } else {
    metric_name
  }
}

# ----------------------------
# helpers: header plots
# ----------------------------
.make_col_header <- function(txt) {
  ggplot() +
    theme_void() +
    annotate("text", x = 0, y = 0, label = txt, hjust = 0.5, size = 6) +
    coord_cartesian(clip = "off")
}

# ----------------------------
# helper: build one ECDF panel (no facet)
# ----------------------------
.build_ecdf_panel <- function(dat, metric_info, lt_vals, show_legend = FALSE) {
  p <- ggplot(dat, aes(x = metric_value, color = embedding, linetype = linetype_by_n)) +
    stat_ecdf(geom = "step", linewidth = 0.9, alpha = 0.95) +
    base_theme +
    theme(
      legend.position = if (show_legend) "right" else "none",
      #axis.text.x  = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18)
    ) +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = lt_vals) +
    labs(color = "Embedding", linetype = "Embedding")

  if (isTRUE(metric_info$is_pvalue)) {
    p <- p + geom_vline(xintercept = 0.05, color = "black", alpha = 0.35, linetype = "dashed")
  }
  p
}
lt_vals <- c("fixed" = "solid", linetype_map)

# ----------------------------
# build a 3-panel row (T1E | Power | Power zoom)
# base_metric_info: if provided, use split column for varying embeddings and base column for fixed embeddings
# ----------------------------
.make_zoom_row <- function(metric_name, metric_info, zoom_max = 0.1, show_legend_in = c("power"), y_label = "ECDF",
                           base_metric_info = NULL) {

  metric_col <- metric_info$col
  
  # Conditional column selection: if base_metric_info provided (split case),
  # use split column for trained embeddings, base column for others
  if (!is.null(base_metric_info)) {
    base_col <- base_metric_info$col
    ecdf_data <- plot_data_diag %>%
      dplyr::mutate(
        metric_value = ifelse(
          as.character(embedding) %in% varying_embeddings_disp,
          !!rlang::sym(metric_col),
          !!rlang::sym(base_col)
        )
      ) %>%
      dplyr::select(seed, n_sample, condition, embedding, metric_value) %>%
      dplyr::filter(!is.na(metric_value)) %>%
      dplyr::filter(
        (as.character(embedding) %in% varying_embeddings_disp) |
          (as.character(embedding) %in% fixed_embeddings_disp &
             as.character(n_sample) == as.character(n_ref_for_fixed))
      ) %>%
      dplyr::mutate(
        linetype_by_n = ifelse(as.character(embedding) %in% varying_embeddings_disp,
                               as.character(n_sample),
                               "fixed")
      )
  } else {
    # Original single-column logic
    ecdf_data <- plot_data_diag %>%
      dplyr::select(seed, n_sample, condition, embedding, metric_value = !!rlang::sym(metric_col)) %>%
      dplyr::filter(!is.na(metric_value)) %>%
      dplyr::filter(
        (as.character(embedding) %in% varying_embeddings_disp) |
          (as.character(embedding) %in% fixed_embeddings_disp &
             as.character(n_sample) == as.character(n_ref_for_fixed))
      ) %>%
      dplyr::mutate(
        linetype_by_n = ifelse(as.character(embedding) %in% varying_embeddings_disp,
                               as.character(n_sample),
                               "fixed")
      )
  }

  # floor p-values (important for clean ECDF + future log options)
  if (isTRUE(metric_info$is_pvalue)) {
    ecdf_data <- ecdf_data %>%
      dplyr::mutate(metric_value = ifelse(metric_value < min_pval_plot | metric_value == 0,
                                          min_pval_plot, metric_value))
  }

  dat_t1e   <- ecdf_data %>% dplyr::filter(condition == "T1E")
  dat_power <- ecdf_data %>% dplyr::filter(condition == "Power")

  # Decide where the legend lives for this row (so patchwork can collect once)
  show_leg_t1e   <- "t1e"   %in% show_legend_in
  show_leg_power <- "power" %in% show_legend_in
  show_leg_zoom  <- "zoom"  %in% show_legend_in

  p_t1e <- .build_ecdf_panel(dat_t1e, metric_info, lt_vals, show_legend = show_leg_t1e) +
    labs(x = NULL, y = y_label) #+
    #theme(
    #  axis.title.y = element_text(size = 14)
    #)

  p_power <- .build_ecdf_panel(dat_power, metric_info, lt_vals, show_legend = show_leg_power) +
    labs(x = NULL, y = NULL)

  p_zoom <- .build_ecdf_panel(dat_power, metric_info, lt_vals, show_legend = show_leg_zoom) +
    labs(x = NULL, y = NULL)
  
  # Apply log10 scale to zoom panel for p-value metrics
  if (isTRUE(metric_info$is_pvalue)) {
    p_zoom <- p_zoom +
      scale_x_log10(
        breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.1),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        limits = c(min_pval_plot, zoom_max)
      )
  } else {
    p_zoom <- p_zoom + coord_cartesian(xlim = c(0, zoom_max))
  }

  # Remove duplicated y-axis labels on middle/right panels
  p_power <- p_power + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  p_zoom  <- p_zoom  + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), aspect.ratio = 1)

  (p_t1e | p_power | p_zoom) + plot_layout(widths = c(1.15, 1.15, 0.9))
}

# ----------------------------
# build combined (pcm on top, rcot below) with shared headers + shared x label
# ----------------------------

# ensure lt_vals exists (same as in your current ECDF code)
# lt_vals <- c("fixed" = "solid", linetype_map)

stopifnot("pcm" %in% names(diagnostic_metrics))
stopifnot("rcot" %in% names(diagnostic_metrics))
stopifnot("pcm_split" %in% names(diagnostic_metrics))
stopifnot("rcot_split" %in% names(diagnostic_metrics))

header_row <- (.make_col_header("T1E") | .make_col_header("Power") | .make_col_header("Power (zoom: x in [0, 0.1])")) +
  plot_layout(widths = c(1.15, 1.15, 0.9), heights=c(1.15, 1.15, 0.9))

# For split metrics: use split column for trained embeddings (Scratch, MedicalNet-ft),
# use non-split column for other embeddings (Freesurfer, cVAE, MedicalNet, FAST)
pcm_row  <- .make_zoom_row("pcm",  diagnostic_metrics[["pcm_split"]],  zoom_max = 0.1, show_legend_in = "power", y_label = "PCM ECDF",
                           base_metric_info = diagnostic_metrics[["pcm"]])
rcot_row <- .make_zoom_row("rcot", diagnostic_metrics[["rcot_split"]], zoom_max = 0.1, show_legend_in = character(0), y_label = "RCoT ECDF",
                           base_metric_info = diagnostic_metrics[["rcot"]])

combined_pcm_rcot_zoom <- (header_row / pcm_row / rcot_row) +
  plot_layout(heights = c(0.12, 1, 1), guides = "collect") &
  theme(legend.position = "right")

# Add ONE shared x-axis label: easiest is to show x title only on bottom row panels
# We'll turn on x-axis title for bottom row only:
combined_pcm_rcot_zoom <- combined_pcm_rcot_zoom &
  theme(axis.title.x = element_blank())

# Rebuild bottom row with x label kept (only once) by editing rcot_row panels:
# (simple approach: add label to the middle plot of rcot row via plot_annotation caption)
combined_pcm_rcot_zoom <- combined_pcm_rcot_zoom +
  plot_annotation(
    caption = "Selection criterion's p-value on validation split",
    theme = theme(plot.caption = element_text(size = 20, hjust = 0.5))
  )

# Save
pdf_path <- file.path(figures_dir, sprintf("diagnostic_ecdf_zoom_pcm_rcot_combined_%d_%d_eps%.1f.pdf", min(seeds), max(seeds), eps_sigmaY))
png_path <- file.path(figures_dir, sprintf("diagnostic_ecdf_zoom_pcm_rcot_combined_%d_%d_eps%.1f.png", min(seeds), max(seeds), eps_sigmaY))

ggsave(pdf_path, combined_pcm_rcot_zoom, width = 18, height = 10, units = "in")
ggsave(png_path, combined_pcm_rcot_zoom, width = 18, height = 10, units = "in", dpi = 300)

message(sprintf("Saved combined PCM+RCoT zoom plot: %s", basename(pdf_path)))


# ============================================================================
# Helper Function: Generate Selection Utility ECDF Plot
# ============================================================================
.generate_selection_utility_ecdf <- function(cit_name, diag_metric_name, 
                                              exclude_candidates,
                                              candidate_embeddings,
                                              alpha = 0.05,
                                              y_label = "ECDF",
                                              verbose = FALSE) {
  
  if (verbose) {
    cat(sprintf("\n--- Generating selection utility for %s CIT with %s diagnostic ---\n", 
                cit_name, diag_metric_name))
  }
  
  # Get diagnostic metric info
  stopifnot(diag_metric_name %in% names(diagnostic_metrics))
  sel_info <- diagnostic_metrics[[diag_metric_name]]
  sel_col  <- sel_info$col
  
  # Determine if using split metric and get base metric
  is_split_metric <- grepl("_split$", diag_metric_name)
  if (is_split_metric) {
    base_metric_name <- .get_base_metric(diag_metric_name)
    base_info <- diagnostic_metrics[[base_metric_name]]
    base_col <- base_info$col
  } else {
    base_col <- sel_col
  }
  
  # Build diagnostic score: higher is better
  diag_sel <- plot_data_diag %>%
    mutate(
      metric_value = ifelse(
        as.character(embedding) %in% varying_embeddings_disp,
        !!sym(sel_col),
        !!sym(base_col)
      )
    ) %>%
    select(seed, n_sample, condition, embedding, metric_value) %>%
    filter(!is.na(metric_value)) %>%
    mutate(
      metric_value = ifelse(isTRUE(sel_info$is_pvalue) & (metric_value < min_pval_plot | metric_value == 0),
                            min_pval_plot, metric_value),
      diag_score = if (isTRUE(sel_info$is_pvalue)) -log10(metric_value) else metric_value
    ) %>%
    filter(embedding %in% candidate_embeddings)
  
  # DNCIT p-values for the chosen CIT
  pvals_sel <- cit_pvals_for_ranks %>%
    filter(cit_label == cit_name) %>%
    filter(embedding %in% candidate_embeddings) %>%
    select(seed, n_sample, condition, embedding, p_value)
  
  # Join CIT p-values with diagnostic scores
  sel_join <- pvals_sel %>%
    left_join(diag_sel %>% select(seed, n_sample, condition, embedding, diag_score),
              by = c("seed", "n_sample", "condition", "embedding")) %>%
    filter(!is.na(diag_score), !is.na(p_value))
  
  if (nrow(sel_join) == 0) {
    warning(sprintf("No data for %s CIT with %s diagnostic", cit_name, diag_metric_name))
    return(NULL)
  }
  
  # Define fixed baseline
  fixed_baseline <- if (cit_name == "RCoT") "MedicalNet" else
    if (cit_name == "PCM") "cVAE" else candidate_embeddings[1]
  
  # Compute selected p-values per strategy
  selection_perf <- sel_join %>%
    group_by(seed, n_sample, condition) %>%
    summarise(
      p_oracle = min(p_value, na.rm = TRUE),
      p_diag = {
        idx <- which.max(diag_score)
        p_value[idx]
      },
      embedding_diag = {
        idx <- which.max(diag_score)
        as.character(embedding[idx])
      },
      p_fixed = p_value[which(embedding == fixed_baseline)[1]],
      p_rand = {
        avail_embeddings <- unique(as.character(embedding))
        if (length(avail_embeddings) > 0) {
          pick <- sample(avail_embeddings, 1)
          p_value[which(as.character(embedding) == pick)[1]]
        } else {
          NA
        }
      },
      .groups = "drop"
    )
  
  # Filter to Power condition only
  selection_perf <- selection_perf %>%
    filter(condition == "Power")
  
  # Convert to long format
  selection_perf <- selection_perf %>%
    as.data.table() %>%
    data.table::melt(
      id.vars = c("seed", "n_sample", "condition"),
      measure.vars = c("p_oracle", "p_diag", "p_fixed", "p_rand"),
      variable.name = "strategy",
      value.name = "p_sel"
    ) %>%
    as.data.frame() %>%
    mutate(
      strategy = recode(strategy,
                        p_oracle = "Oracle",
                        p_diag   = "Selection criterion",
                        p_fixed  = paste0("Fixed (", fixed_baseline, ")"),
                        p_rand   = "Random"),
      minuslog10p = -log10(pmax(p_sel, min_pval_plot)),
      reject = as.integer(p_sel <= alpha)
    )
  
  # ECDF version: Power condition only with log scale
  selection_perf_ecdf <- selection_perf %>%
    mutate(
      p_sel_plot = pmax(p_sel, min_pval_plot)
    )
  
  # Get strategy colors
  strategy_order <- c("Oracle", "Selection criterion", paste0("Fixed (", fixed_baseline, ")"), "Random")
  strategy_colors_vec <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")[1:4]
  strategy_colors <- setNames(strategy_colors_vec, strategy_order)
  
  # Create ECDF plot
  p_ecdf <- ggplot(selection_perf_ecdf, aes(x = p_sel_plot, color = strategy)) +
    stat_ecdf(geom = "step", linewidth = 0.9, alpha = 0.95) +
    facet_wrap(~ n_sample, nrow = 1) +
    scale_x_log10(
      breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3,  0.1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(min_pval_plot, 1)
    ) +
    geom_vline(xintercept = 0.05, color = "black", alpha = 0.4, linetype = "dashed") +
    labs(
      x = NULL,  # Will be added to combined plot
      y = y_label,
      color = "Strategy"
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "right"
    ) +
    scale_color_manual(values = strategy_colors)
  
  return(p_ecdf)
}

# ============================================================================
# NEW FIGURE 2: Utility of diagnostic-based embedding selection
#   - Exclude FAST and Freesurfer from candidate set
#   - Compare Oracle vs Diagnostic-selected vs Fixed vs Random
# ============================================================================
cat("\n=== Creating Selection Utility Plots ===\n")

set.seed(1)

alpha <- 0.05
exclude_candidates <- c("FAST", "Freesurfer")#c("FAST")
candidate_embeddings <- setdiff(levels(cit_pvals_for_ranks$embedding), exclude_candidates)

# Choose which DNCIT CIT to evaluate in Figure 2
# Options present in your script: "PCM" and "RCoT"
cit_for_selection <- "RCoT"  # <-- change to "RCoT" if desired

# Choose which diagnostic drives the selection
diag_metric_for_selection <- "rcot_split"  # <-- e.g. "ftest", "pcm", "r2_residual", ...

stopifnot(diag_metric_for_selection %in% names(diagnostic_metrics))
sel_info <- diagnostic_metrics[[diag_metric_for_selection]]
sel_col  <- sel_info$col

# Determine if using split metric and get base metric
# For split metrics: use split column for trained embeddings (Scratch, MedicalNet-ft),
# use non-split column for other embeddings (Freesurfer, cVAE, MedicalNet, FAST)
is_split_metric <- grepl("_split$", diag_metric_for_selection)
if (is_split_metric) {
  base_metric_name <- .get_base_metric(diag_metric_for_selection)
  base_info <- diagnostic_metrics[[base_metric_name]]
  base_col <- base_info$col
} else {
  base_col <- sel_col
}

# Build diagnostic score: higher is better
# Conditional column selection based on embedding type
diag_sel <- plot_data_diag %>%
  mutate(
    metric_value = ifelse(
      as.character(embedding) %in% varying_embeddings_disp,
      !!sym(sel_col),
      !!sym(base_col)
    )
  ) %>%
  select(seed, n_sample, condition, embedding, metric_value) %>%
  filter(!is.na(metric_value)) %>%
  mutate(
    metric_value = ifelse(isTRUE(sel_info$is_pvalue) & (metric_value < min_pval_plot | metric_value == 0),
                          min_pval_plot, metric_value),
    diag_score = if (isTRUE(sel_info$is_pvalue)) -log10(metric_value) else metric_value
  ) %>%
  filter(embedding %in% candidate_embeddings)

# DNCIT p-values for the chosen CIT (PCM or RCoT) and candidates
pvals_sel <- cit_pvals_for_ranks %>%
  filter(cit_label == cit_for_selection) %>%
  filter(embedding %in% candidate_embeddings) %>%
  select(seed, n_sample, condition, embedding, p_value)

cat("\n=== Selection Data Diagnostics ===\n")
cat(sprintf("CIT: %s, Diagnostic metric: %s\n", cit_for_selection, diag_metric_for_selection))
cat(sprintf("\nEmbeddings with CIT p-values: %s\n", 
            paste(sort(unique(as.character(pvals_sel$embedding))), collapse = ", ")))
cat(sprintf("Count by embedding:\n"))
print(table(pvals_sel$embedding))

cat(sprintf("\nEmbeddings with diagnostic scores: %s\n", 
            paste(sort(unique(as.character(diag_sel$embedding))), collapse = ", ")))
cat(sprintf("Count by embedding:\n"))
print(table(diag_sel$embedding))

cat(sprintf("\nRows in pvals_sel: %d\n", nrow(pvals_sel)))
cat(sprintf("Rows in diag_sel: %d\n", nrow(diag_sel)))

# Join CIT p-values with diagnostic scores
# NOTE: This join requires diagnostic scores for all embeddings. If an embedding
# (e.g., Scratch) has CIT p-values but no diagnostic scores, it will be EXCLUDED
# from ALL selection strategies (Oracle, Diagnostic, Fixed, Random).
# If diagnostic output shows missing embeddings, either:
#   1. Compute diagnostics for those embeddings, OR
#   2. Exclude them from candidate_embeddings
sel_join <- pvals_sel %>%
  left_join(diag_sel %>% select(seed, n_sample, condition, embedding, diag_score),
            by = c("seed", "n_sample", "condition", "embedding")) %>%
  filter(!is.na(diag_score), !is.na(p_value))

cat(sprintf("\nRows after join and filter: %d\n", nrow(sel_join)))
cat(sprintf("Embeddings in joined data: %s\n", 
            paste(sort(unique(as.character(sel_join$embedding))), collapse = ", ")))
if (nrow(sel_join) > 0) {
  cat(sprintf("Count by embedding:\n"))
  print(table(sel_join$embedding))
}

# Check for embeddings with p-values but no diagnostic scores
missing_diag_embeddings <- setdiff(
  unique(as.character(pvals_sel$embedding)),
  unique(as.character(sel_join$embedding))
)
if (length(missing_diag_embeddings) > 0) {
  cat(sprintf("\nWARNING: The following embeddings have CIT p-values but are MISSING diagnostic scores:\n"))
  cat(sprintf("  %s\n", paste(missing_diag_embeddings, collapse = ", ")))
  cat(sprintf("These embeddings are EXCLUDED from selection analysis (including Oracle).\n"))
  cat(sprintf("To include them, compute diagnostic scores for these embeddings.\n"))
}
cat("===\n\n")

if (nrow(sel_join) == 0) {
  warning("Figure 2 selection join produced no rows (check metric choice / missing diagnostics).")
} else {

  # Define a fixed baseline (prefer cVAE, else MedicalNet, else first available)
  fixed_baseline <- if ("RCoT" %in% cit_for_selection) "MedicalNet" else
    if ("PCM" %in% cit_for_selection) "cVAE" else candidate_embeddings[1]

  # Compute selected p-values per strategy for each (seed, n_sample, condition)
  # Also track which embedding is selected by diagnostic strategy
  selection_perf <- sel_join %>%
    group_by(seed, n_sample, condition) %>%
    summarise(
      p_oracle = min(p_value, na.rm = TRUE),
      p_diag = {
        idx <- which.max(diag_score)
        p_value[idx]
      },
      embedding_diag = {
        idx <- which.max(diag_score)
        as.character(embedding[idx])
      },
      p_fixed = p_value[which(embedding == fixed_baseline)[1]],
      p_rand = {
        # Fix: sample from actual embedding values in this group, not character conversion
        avail_embeddings <- unique(as.character(embedding))
        if (length(avail_embeddings) > 0) {
          pick <- sample(avail_embeddings, 1)
          p_value[which(as.character(embedding) == pick)[1]]
        } else {
          NA
        }
      },
      .groups = "drop"
    )
  
  # Count how often each embedding is selected by diagnostic strategy
  diag_selection_counts <- selection_perf %>%
    count(embedding_diag, n_sample, condition, name = "count") %>%
    group_by(n_sample, condition) %>%
    mutate(
      total = sum(count),
      proportion = count / total
    ) %>%
    ungroup() %>%
    arrange(condition, n_sample, desc(count))
  
  cat("\n=== Diagnostic Embedding Selection Counts ===\n")
  cat(sprintf("Diagnostic metric: %s\n", diag_metric_for_selection))
  cat(sprintf("CIT: %s\n", cit_for_selection))
  cat("\nSelection frequencies:\n")
  print(diag_selection_counts)
  
  # Overall counts across all conditions and sample sizes
  diag_selection_overall <- selection_perf %>%
    count(embedding_diag, name = "count") %>%
    mutate(
      total = sum(count),
      proportion = count / total
    ) %>%
    arrange(desc(count))
  
  cat("\nOverall selection frequencies (across all conditions and sample sizes):\n")
  print(diag_selection_overall)
  
  # Save to CSV
  output_file_selection <- file.path(figures_dir, sprintf("diagnostic_selection_counts_%s_by_%s_%d_%d.csv",
                                                           cit_for_selection, diag_metric_for_selection, min(seeds), max(seeds)))
  fwrite(diag_selection_counts, output_file_selection)
  cat(sprintf("\nSaved selection counts to: %s\n", basename(output_file_selection)))
  
  # Analyze Power condition: when Freesurfer is not selected, how often was the selected embedding optimal?
  power_analysis <- selection_perf %>%
    filter(condition == "Power") %>%
    # Join with sel_join to get all embeddings and their p-values for each (seed, n_sample)
    left_join(
      sel_join %>% 
        filter(condition == "Power") %>%
        group_by(seed, n_sample) %>%
        summarise(
          min_p_value = min(p_value, na.rm = TRUE),
          embedding_with_min_p = embedding[which.min(p_value)][1],
          .groups = "drop"
        ),
      by = c("seed", "n_sample")
    ) %>%
    mutate(
      freesurfer_not_selected = (embedding_diag != "Freesurfer"),
      selected_was_optimal = (embedding_diag == as.character(embedding_with_min_p))
    )
  
  # Count how often Freesurfer is not selected
  freesurfer_not_selected_count <- power_analysis %>%
    summarise(
      total_cases = n(),
      freesurfer_not_selected = sum(freesurfer_not_selected, na.rm = TRUE),
      freesurfer_not_selected_prop = mean(freesurfer_not_selected, na.rm = TRUE)
    )
  
  # When Freesurfer is not selected, count how often the selected embedding had the minimum p-value
  optimal_when_not_freesurfer <- power_analysis %>%
    filter(freesurfer_not_selected) %>%
    summarise(
      total_when_not_freesurfer = n(),
      selected_was_optimal_count = sum(selected_was_optimal, na.rm = TRUE),
      selected_was_optimal_prop = mean(selected_was_optimal, na.rm = TRUE)
    )
  
  cat("\n=== Power Condition Analysis: Freesurfer Non-Selection ===\n")
  cat(sprintf("Total cases in Power condition: %d\n", freesurfer_not_selected_count$total_cases))
  cat(sprintf("Freesurfer NOT selected: %d (%.2f%%)\n", 
              freesurfer_not_selected_count$freesurfer_not_selected,
              freesurfer_not_selected_count$freesurfer_not_selected_prop * 100))
  cat(sprintf("\nWhen Freesurfer was NOT selected:\n"))
  cat(sprintf("  Total cases: %d\n", optimal_when_not_freesurfer$total_when_not_freesurfer))
  cat(sprintf("  Selected embedding had minimum p-value: %d (%.2f%%)\n",
              optimal_when_not_freesurfer$selected_was_optimal_count,
              optimal_when_not_freesurfer$selected_was_optimal_prop * 100))
  
  # Detailed breakdown by sample size
  optimal_by_sample <- power_analysis %>%
    filter(freesurfer_not_selected) %>%
    group_by(n_sample) %>%
    summarise(
      total = n(),
      selected_was_optimal = sum(selected_was_optimal, na.rm = TRUE),
      prop_optimal = mean(selected_was_optimal, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\nBreakdown by sample size (when Freesurfer not selected):\n")
  print(optimal_by_sample)
  
  # Save detailed results
  # Include the p-value of the selected embedding (p_diag is already in power_analysis)
  power_analysis_detailed <- power_analysis %>%
    select(seed, n_sample, embedding_diag, p_diag, min_p_value, embedding_with_min_p, 
           freesurfer_not_selected, selected_was_optimal)
  
  output_file_power <- file.path(figures_dir, sprintf("power_freesurfer_analysis_%s_by_%s_%d_%d.csv",
                                                       cit_for_selection, diag_metric_for_selection, min(seeds), max(seeds)))
  fwrite(power_analysis_detailed, output_file_power)
  cat(sprintf("\nSaved detailed Power analysis to: %s\n", basename(output_file_power)))
  
  # Summary of oracle p-values by sample size for Power condition (before filtering)
  oracle_summary <- selection_perf %>%
    filter(condition == "Power") %>%
    group_by(n_sample) %>%
    summarise(
      oracle_min = min(p_oracle, na.rm = TRUE),
      oracle_q25 = quantile(p_oracle, 0.25, na.rm = TRUE),
      oracle_median = median(p_oracle, na.rm = TRUE),
      oracle_q75 = quantile(p_oracle, 0.75, na.rm = TRUE),
      oracle_max = max(p_oracle, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n=== Oracle P-value Summary (Power) ===\n")
  print(oracle_summary)
  cat("\n")
  
  # Filter to Power condition only
  selection_perf <- selection_perf %>%
    filter(condition == "Power")
  
  # Convert to long format
  selection_perf <- selection_perf %>%
    as.data.table() %>%
    data.table::melt(
      id.vars = c("seed", "n_sample", "condition"),
      measure.vars = c("p_oracle", "p_diag", "p_fixed", "p_rand"),
      variable.name = "strategy",
      value.name = "p_sel"
    ) %>%
    as.data.frame() %>%
    mutate(
      strategy = recode(strategy,
                        p_oracle = "Oracle",
                        p_diag   = "Diagnostic-selected",
                        p_fixed  = paste0("Fixed (", fixed_baseline, ")"),
                        p_rand   = "Random"),
      minuslog10p = -log10(pmax(p_sel, min_pval_plot)),
      reject = as.integer(p_sel <= alpha)
    )

  # Boxplots: -log10(p) for Power condition
  p_fig2 <- ggplot(selection_perf, aes(x = strategy, y = minuslog10p)) +
    geom_boxplot(outlier.alpha = 0.15) +
    facet_wrap(~ n_sample, nrow = 1) +
    labs(
      x = NULL,
      y = "-log10(p-value)",
      title = sprintf("Embedding selection utility - Power (%s CIT, %s diagnostic)", cit_for_selection, diag_metric_for_selection),
      subtitle = sprintf("Candidates exclude %s", paste(exclude_candidates, collapse = " & "))
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "none"
    )
  
  # Add y-axis labels manually via annotation
  # Extract the plot to add custom y-axis labels
  # T1E panels should show "p-value" and Power panels should show "-log10(p)"
  # This is handled by the subtitle, but we can also add it to strip labels

  pdf_path2 <- file.path(figures_dir, sprintf("selection_utility_%s_by_%s_%d_%d_eps%.1f_exclude_%s.pdf",
                                              cit_for_selection, diag_metric_for_selection, min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))
  png_path2 <- file.path(figures_dir, sprintf("selection_utility_%s_by_%s_%d_%d_eps%.1f_exclude_%s.png",
                                              cit_for_selection, diag_metric_for_selection, min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))

  ggsave(pdf_path2, p_fig2, width = 14, height = 7, units = "in")
  ggsave(png_path2, p_fig2, width = 14, height = 7, units = "in", dpi = 300)

  message(sprintf("Saved selection utility plot: %s", basename(pdf_path2)))
  
  # ECDF version: Power condition only with log scale
  # Floor p-values for log plotting
  selection_perf_ecdf <- selection_perf %>%
    mutate(
      p_sel_plot = pmax(p_sel, min_pval_plot)
    )
  
  # Get strategy colors from paletteer (using first 4 colors from Classic_10_Medium)
  strategy_order <- c("Oracle", "Diagnostic-selected", paste0("Fixed (", fixed_baseline, ")"), "Random")
  strategy_colors_vec <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")[1:4]
  strategy_colors <- setNames(strategy_colors_vec, strategy_order)
  
  p_fig2_ecdf <- ggplot(selection_perf_ecdf, aes(x = p_sel_plot, color = strategy)) +
    stat_ecdf(geom = "step", linewidth = 0.9, alpha = 0.95) +
    facet_wrap(~ n_sample, nrow = 1) +
    scale_x_log10(
    	  breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3,  0.1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(min_pval_plot, 1)
    ) +
    geom_vline(xintercept = 0.05, color = "black", alpha = 0.4, linetype = "dashed") +
    labs(
      x = "p-value",
      y = "ECDF",
      color = "Strategy",
      title = sprintf("Embedding selection utility - Power (%s CIT, %s diagnostic)", cit_for_selection, diag_metric_for_selection),
      subtitle = sprintf("Candidates exclude %s", paste(exclude_candidates, collapse = " & "))
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "right"
    ) +
    scale_color_manual(values = strategy_colors)
  
  pdf_path2_ecdf <- file.path(figures_dir, sprintf("selection_utility_ecdf_%s_by_%s_%d_%d_eps%.1f_exclude_%s.pdf",
                                              cit_for_selection, diag_metric_for_selection, min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))
  png_path2_ecdf <- file.path(figures_dir, sprintf("selection_utility_ecdf_%s_by_%s_%d_%d_eps%.1f_exclude_%s.png",
                                              cit_for_selection, diag_metric_for_selection, min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))

  ggsave(pdf_path2_ecdf, p_fig2_ecdf, width = 14, height = 7, units = "in")
  ggsave(png_path2_ecdf, p_fig2_ecdf, width = 14, height = 7, units = "in", dpi = 300)

  message(sprintf("Saved selection utility ECDF plot: %s", basename(pdf_path2_ecdf)))
}

cat("\n=== Selection Utility Plots Complete ===\n")

# ============================================================================
# Generate Combined Selection Utility Plots (PCM + RCoT)
# ============================================================================
cat("\n=== Creating Combined Selection Utility Plots ===\n")

set.seed(1)
alpha <- 0.05
exclude_candidates <- c('FAST')#c("FAST", "Freesurfer")
candidate_embeddings <- setdiff(levels(cit_pvals_for_ranks$embedding), exclude_candidates)

# Generate individual plots
cat("\n--- Generating individual ECDF plots ---\n")
plot_pcm_regular <- .generate_selection_utility_ecdf("PCM", "pcm", exclude_candidates, candidate_embeddings, 
                                                      alpha = alpha, y_label = "ECDF (PCM)", verbose = TRUE)
plot_rcot_regular <- .generate_selection_utility_ecdf("RCoT", "rcot", exclude_candidates, candidate_embeddings, 
                                                       alpha = alpha, y_label = "ECDF (RCoT)", verbose = TRUE)
plot_pcm_split <- .generate_selection_utility_ecdf("PCM", "pcm_split", exclude_candidates, candidate_embeddings, 
                                                    alpha = alpha, y_label = "ECDF (PCM)", verbose = TRUE)
plot_rcot_split <- .generate_selection_utility_ecdf("RCoT", "rcot_split", exclude_candidates, candidate_embeddings, 
                                                     alpha = alpha, y_label = "ECDF (RCoT)", verbose = TRUE)

# Combine plots using patchwork
cat("\n--- Combining plots ---\n")

# Figure 1: Regular metrics (PCM + RCoT)
if (!is.null(plot_pcm_regular) && !is.null(plot_rcot_regular)) {
  # Remove x-axis elements from top plot, add x-axis label to bottom plot
  plot_pcm_regular_mod <- plot_pcm_regular + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  plot_rcot_regular_mod <- plot_rcot_regular + 
    labs(x = "p-value")
  
  combined_regular <- plot_pcm_regular_mod / plot_rcot_regular_mod +
    plot_layout(guides = "collect", axes = "collect_x") &
    theme(legend.position = "right")
  
  # Save combined regular
  pdf_path_regular <- file.path(figures_dir, sprintf("selection_utility_ecdf_combined_regular_%d_%d_eps%.1f_exclude_%s.pdf",
                                                      min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))
  png_path_regular <- file.path(figures_dir, sprintf("selection_utility_ecdf_combined_regular_%d_%d_eps%.1f_exclude_%s.png",
                                                      min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))
  
  ggsave(pdf_path_regular, combined_regular, width = 14, height = 12, units = "in")
  ggsave(png_path_regular, combined_regular, width = 14, height = 12, units = "in", dpi = 300)
  
  cat(sprintf("Saved combined regular plot: %s\n", basename(pdf_path_regular)))
}

# Figure 2: Split metrics (PCM + RCoT)
if (!is.null(plot_pcm_split) && !is.null(plot_rcot_split)) {
  # Remove x-axis elements from top plot, add x-axis label to bottom plot
  plot_pcm_split_mod <- plot_pcm_split + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  plot_rcot_split_mod <- plot_rcot_split + 
    labs(x = "DNCIT's p-value on test split")
  
  combined_split <- plot_pcm_split_mod / plot_rcot_split_mod +
    plot_layout(guides = "collect", axes = "collect_x") &
    theme(legend.position = "right")
  
  # Save combined split
  pdf_path_split <- file.path(figures_dir, sprintf("selection_utility_ecdf_combined_split_%d_%d_eps%.1f_exclude_%s.pdf",
                                                    min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))
  png_path_split <- file.path(figures_dir, sprintf("selection_utility_ecdf_combined_split_%d_%d_eps%.1f_exclude_%s.png",
                                                    min(seeds), max(seeds), eps_sigmaY, paste(exclude_candidates, collapse = "_")))
  
  ggsave(pdf_path_split, combined_split, width = 14, height = 12, units = "in")
  ggsave(png_path_split, combined_split, width = 14, height = 12, units = "in", dpi = 300)
  
  cat(sprintf("Saved combined split plot: %s\n", basename(pdf_path_split)))
}

cat("\n=== Combined Selection Utility Plots Complete ===\n")

cat("\n=== Complete ===\n")



