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
n_samples <- c(460, 1100, 5000)#, 10000)
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
    
    sample_size_cols <- names(pval_df)[1:min(length(n_samples), ncol(pval_df))]
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
                          "f_test_pvalue", "globaltest_pvalue", "dcor_pvalue")
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
         r2_test, f_test_pvalue, globaltest_pvalue, dcor_pvalue) %>%
  as.data.table()

# Remove duplicates by taking the first row for each unique combination
diagnostics_summary <- diagnostics_summary[, .SD[1], 
                                          by = .(seed, n_sample, condition, embedding_name, y_type)]

# Reshape to have separate columns for original and residual R^2
r2_original <- diagnostics_summary[y_type == "original", 
                                   .(seed, n_sample, condition, embedding_name, r2_original = r2_test)]
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

cat("\n=== Complete ===\n")

