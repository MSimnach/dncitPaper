### Correlation between CIT P-values and Diagnostics
# Compute correlations between CIT p-values and diagnostic metrics
# (R^2 original, R^2 residual, F-test, globaltest, dcor.test p-values)

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Configuration
seeds <- c(601:625)
n_samples <- c(460, 1100, 5000, 10000)
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
  pval_dir <- file.path(results_base_path, cond, "p-values", "seeds_601_625")
  
  if (!dir.exists(pval_dir)) {
    warning(sprintf("Directory not found: %s", pval_dir))
    next
  }
  
  csv_files <- list.files(pval_dir, pattern = "\\.csv$", full.names = TRUE)
  cat(sprintf("Condition: %s - Found %d CSV files\n", cond, length(csv_files)))
  
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
        cor_results <- data.frame(
          cit = cit_name,
          embedding = emb,
          condition = cond,
          n_sample = n_samp,
          cor_r2_orig = cor(merged_data$p_value, merged_data$r2_original, use = "complete.obs"),
          cor_r2_resid = cor(merged_data$p_value, merged_data$r2_residual, use = "complete.obs"),
          cor_ftest = cor(merged_data$p_value, merged_data$f_test_pvalue, use = "complete.obs"),
          cor_globaltest = cor(merged_data$p_value, merged_data$globaltest_pvalue, use = "complete.obs"),
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
      cor_results <- data.frame(
        cit = cit_name,
        condition = cond,
        n_sample = n_samp,
        cor_r2_orig = cor(merged_data$p_value, merged_data$r2_original, use = "complete.obs"),
        cor_r2_resid = cor(merged_data$p_value, merged_data$r2_residual, use = "complete.obs"),
        cor_ftest = cor(merged_data$p_value, merged_data$f_test_pvalue, use = "complete.obs"),
        cor_globaltest = cor(merged_data$p_value, merged_data$globaltest_pvalue, use = "complete.obs"),
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
# Save Results
# ============================================================================
cat("=== Saving Results ===\n")

output_dir <- "/sc/home/marco.simnacher/dncitPaper/inst/diagnostic_embedding"

# Save per-embedding correlations
output_file_per_emb <- file.path(output_dir, "correlation_cit_diagnostics_per_embedding.csv")
fwrite(correlation_df, output_file_per_emb)
cat(sprintf("Saved per-embedding correlations: %s\n", output_file_per_emb))

# Save across-all-embeddings correlations
output_file_pooled <- file.path(output_dir, "correlation_cit_diagnostics_pooled.csv")
fwrite(correlation_pooled_df, output_file_pooled)
cat(sprintf("Saved pooled correlations: %s\n", output_file_pooled))

cat("\n=== Summary ===\n")
cat("Per-embedding correlations:\n")
print(head(correlation_df))
cat("\nPooled correlations:\n")
print(head(correlation_pooled_df))

# ============================================================================
# Create Diagnostic Boxplots
# ============================================================================
cat("\n=== Creating Diagnostic Boxplots ===\n")

# Define embeddings
constant_embeddings <- c("fastsurfer", "freesurfer", "condVAE", "medicalnet")
varying_embeddings <- c('medicalnet_ft', 'scratch')
available_embeddings <- unique(combined_pvals$embedding)

constant_embeddings_plot <- intersect(constant_embeddings, available_embeddings)
varying_embeddings_plot <- intersect(varying_embeddings, available_embeddings)
all_embeddings_plot <- c(constant_embeddings_plot, varying_embeddings_plot)

# Create color palette
n_colors <- length(all_embeddings_plot)
base_colors <- RColorBrewer::brewer.pal(min(8, n_colors), "Set2")
color_palette <- setNames(
  colorRampPalette(base_colors)(n_colors),
  all_embeddings_plot
)

# Merge CIT p-values with diagnostics for plotting
plot_data <- merge(combined_pvals, diagnostics_wide,
                   by.x = c("seed", "n_sample", "condition", "embedding"),
                   by.y = c("seed", "n_sample", "condition", "embedding_name"),
                   all = FALSE)

# Convert to factors for plotting
plot_data <- plot_data %>%
  mutate(
    n_sample = factor(n_sample, levels = n_samples),
    condition = factor(condition, levels = conditions),
    embedding = factor(embedding, levels = all_embeddings_plot),
    cit = factor(cit)
  ) %>%
  filter(embedding %in% all_embeddings_plot)

# Create figures directory
figures_dir <- file.path(output_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum p-value for plotting
min_pval_plot <- 1e-16

# Get unique CITs
available_cits <- unique(plot_data$cit)

# Create plots for each CIT
for (cit_name in available_cits) {
  cat(sprintf("\nCreating diagnostic plots for CIT: %s\n", cit_name))
  
  # Filter data for this CIT
  plot_data_cit <- plot_data %>% filter(cit == cit_name)
  
  if (nrow(plot_data_cit) == 0) {
    warning(sprintf("No data available for CIT: %s", cit_name))
    next
  }
  
  # ===== R^2 Original =====
  p_r2_orig <- plot_data_cit %>%
    ggplot(aes(x = n_sample, y = r2_original, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = expression(R^2 ~ "(Original)"),
      fill = "Embedding",
      title = bquote("R"^2 ~ "Original for" ~ .(cit_name) ~ "(" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold")
    ) +
    scale_fill_manual(values = color_palette) +
    ylim(0, NA)
  
  # Save R^2 Original
  png_path <- file.path(figures_dir, sprintf("diagnostic_r2_original_%s_%d_%d.png", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(png_path, plot = p_r2_orig, width = 16, height = 7, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(png_path)))
  
  pdf_path <- file.path(figures_dir, sprintf("diagnostic_r2_original_%s_%d_%d.pdf", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(pdf_path, plot = p_r2_orig, width = 16, height = 7, units = "in")
  
  # ===== R^2 Residual =====
  p_r2_resid <- plot_data_cit %>%
    ggplot(aes(x = n_sample, y = r2_residual, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = expression(R^2 ~ "(Residual)"),
      fill = "Embedding",
      title = bquote("R"^2 ~ "Residual for" ~ .(cit_name) ~ "(" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold")
    ) +
    scale_fill_manual(values = color_palette) +
    ylim(0, NA)
  
  # Save R^2 Residual
  png_path <- file.path(figures_dir, sprintf("diagnostic_r2_residual_%s_%d_%d.png", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(png_path, plot = p_r2_resid, width = 16, height = 7, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(png_path)))
  
  pdf_path <- file.path(figures_dir, sprintf("diagnostic_r2_residual_%s_%d_%d.pdf", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(pdf_path, plot = p_r2_resid, width = 16, height = 7, units = "in")
  
  # ===== F-test P-value =====
  plot_data_ftest <- plot_data_cit %>%
    mutate(f_test_pvalue = ifelse(f_test_pvalue < min_pval_plot | f_test_pvalue == 0, 
                                   min_pval_plot, f_test_pvalue))
  
  p_ftest <- plot_data_ftest %>%
    ggplot(aes(x = n_sample, y = f_test_pvalue, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = "F-test P-value",
      fill = "Embedding",
      title = bquote("F-test P-values for" ~ .(cit_name) ~ "(" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold")
    ) +
    scale_fill_manual(values = color_palette) +
    scale_y_log10(
      breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.5, 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(min_pval_plot, 1)
    )
  
  # Save F-test
  png_path <- file.path(figures_dir, sprintf("diagnostic_ftest_%s_%d_%d.png", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(png_path, plot = p_ftest, width = 16, height = 7, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(png_path)))
  
  pdf_path <- file.path(figures_dir, sprintf("diagnostic_ftest_%s_%d_%d.pdf", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(pdf_path, plot = p_ftest, width = 16, height = 7, units = "in")
  
  # ===== Globaltest P-value =====
  plot_data_globaltest <- plot_data_cit %>%
    mutate(globaltest_pvalue = ifelse(globaltest_pvalue < min_pval_plot | globaltest_pvalue == 0, 
                                       min_pval_plot, globaltest_pvalue))
  
  p_globaltest <- plot_data_globaltest %>%
    ggplot(aes(x = n_sample, y = globaltest_pvalue, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = "Globaltest P-value",
      fill = "Embedding",
      title = bquote("Globaltest P-values for" ~ .(cit_name) ~ "(" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold")
    ) +
    scale_fill_manual(values = color_palette) +
    scale_y_log10(
      breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.5, 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(min_pval_plot, 1)
    )
  
  # Save Globaltest
  png_path <- file.path(figures_dir, sprintf("diagnostic_globaltest_%s_%d_%d.png", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(png_path, plot = p_globaltest, width = 16, height = 7, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(png_path)))
  
  pdf_path <- file.path(figures_dir, sprintf("diagnostic_globaltest_%s_%d_%d.pdf", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(pdf_path, plot = p_globaltest, width = 16, height = 7, units = "in")
  
  # ===== Dcor P-value =====
  plot_data_dcor <- plot_data_cit %>%
    mutate(dcor_pvalue = ifelse(dcor_pvalue < min_pval_plot | dcor_pvalue == 0, 
                                 min_pval_plot, dcor_pvalue))
  
  p_dcor <- plot_data_dcor %>%
    ggplot(aes(x = n_sample, y = dcor_pvalue, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = "Dcor P-value",
      fill = "Embedding",
      title = bquote("Dcor P-values for" ~ .(cit_name) ~ "(" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(face = "bold")
    ) +
    scale_fill_manual(values = color_palette) +
    scale_y_log10(
      breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.5, 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(min_pval_plot, 1)
    )
  
  # Save Dcor
  png_path <- file.path(figures_dir, sprintf("diagnostic_dcor_%s_%d_%d.png", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(png_path, plot = p_dcor, width = 16, height = 7, units = "in", dpi = 300)
  cat(sprintf("  Saved: %s\n", basename(png_path)))
  
  pdf_path <- file.path(figures_dir, sprintf("diagnostic_dcor_%s_%d_%d.pdf", 
                                             cit_name, min(seeds), max(seeds)))
  ggsave(pdf_path, plot = p_dcor, width = 16, height = 7, units = "in")
}

cat("\n=== Diagnostic Visualization Complete ===\n")
cat(sprintf("All plots saved to: %s\n", figures_dir))

cat("\n=== Complete ===\n")

