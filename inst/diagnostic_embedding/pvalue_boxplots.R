### P-value Extraction and Boxplot Visualization
# Extract CIT p-values from Results/ directory and create boxplots
# showing p-values by CIT, embedding, condition, and sample size

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# ============================================================================
# Configuration Variables
# ============================================================================

# Seeds to analyze
seeds <- c(561:585)

# Sample sizes to include
n_samples <- c(256, 460, 825, 1100, 5000)#c(460, 1100, 5000, 10000)

# Conditions to analyze
conditions <- c("CI", "No_CI")

# Epsilon sigma Y parameter
eps_sigmaY <- 0.5

# Define p-values base directories
results_base_path <- "/sc/home/marco.simnacher/dncitPaper/Results"

# Output directory for plots
output_dir <- "/sc/home/marco.simnacher/dncitPaper/inst/diagnostic_embedding/figures"

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Minimum p-value for plotting (to avoid log10(0) = -Inf)
min_pval_plot <- 1e-16

# Define embedding categories
constant_embeddings <- c("fastsurfer", "freesurfer", "condVAE", "medicalnet")
varying_embeddings <- c('medicalnet_ft', 'scratch')

cat("=== P-value Extraction and Visualization ===\n")
cat(sprintf("Seeds: %d to %d (%d seeds)\n", min(seeds), max(seeds), length(seeds)))
cat(sprintf("Sample sizes: %s\n", paste(n_samples, collapse = ", ")))
cat(sprintf("Conditions: %s\n", paste(conditions, collapse = ", ")))
cat(sprintf("Epsilon sigma Y: %.2f\n\n", eps_sigmaY))

# ============================================================================
# Helper Functions
# ============================================================================

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

# Function to extract CIT name from filename
extract_cit_from_filename <- function(filename) {
  # Pattern: 1_0_[eps]_0_fastsurfer_[embedding]_ukb_z1_squared_[CIT].csv
  # Extract everything between "squared_" and ".csv"
  pattern <- "squared_(.+?)\\.csv"
  match <- regmatches(filename, regexec(pattern, filename))
  if (length(match[[1]]) > 1) {
    return(match[[1]][2])
  } else {
    return(NA)
  }
}

# ============================================================================
# Load P-values from CSV Files
# ============================================================================

cat("\n=== Loading P-values from CIT Tests ===\n")

# Construct seed directory name
seed_dir_name <- sprintf("seeds_%d_%d", min(seeds), max(seeds))

# Initialize list to store p-value data
pval_data_list <- list()

for (cond in conditions) {
  pval_dir <- file.path(results_base_path, cond, "p-values", seed_dir_name)
  
  # Check if directory exists
  if (!dir.exists(pval_dir)) {
    warning(sprintf("Directory not found: %s", pval_dir))
    next
  }
  
  # Get all CSV files in the directory that match the specified eps_sigmaY
  # Pattern: 1_0_[eps_sigmaY]_0_...csv
  pattern_with_eps <- sprintf("^1_0_%d_0_.*\\.csv$", eps_sigmaY)
  csv_files <- list.files(pval_dir, pattern = pattern_with_eps, full.names = TRUE)
  
  cat(sprintf("\nCondition: %s - Found %d CSV files\n", cond, length(csv_files)))
  cat(sprintf("  Directory: %s\n", pval_dir))
  
  for (csv_file in csv_files) {
    filename <- basename(csv_file)
    embedding_name <- extract_embedding_from_filename(filename)
    cit_name <- extract_cit_from_filename(filename)
    
    if (is.na(embedding_name)) {
      warning(sprintf("Could not extract embedding name from: %s", filename))
      next
    }
    
    if (is.na(cit_name)) {
      warning(sprintf("Could not extract CIT name from: %s", filename))
      next
    }
    
    cat(sprintf("  Loading: %s (embedding: %s, CIT: %s)\n", filename, embedding_name, cit_name))
    
    # Read the CSV
    pval_df <- fread(csv_file, nThread = 1)
    
    # Drop the first column by position (contains just 1:25 row indices)
    pval_df <- pval_df[, -1]
    
    # Get column names for sample sizes (should be V1, V2, V3, V4)
    sample_size_cols <- names(pval_df)[1:min(length(n_samples), ncol(pval_df))]
    
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
    pval_long$cit <- cit_name
    
    # Store in list
    pval_data_list[[length(pval_data_list) + 1]] <- pval_long[, .(seed, n_sample, embedding, condition, cit, p_value)]
  }
}

# Combine all p-value data
if (length(pval_data_list) > 0) {
  combined_pvals <- rbindlist(pval_data_list)
  
  cat(sprintf("\n=== Data Summary ===\n"))
  cat(sprintf("Total p-value records loaded: %d\n", nrow(combined_pvals)))
  cat(sprintf("Unique embeddings: %d\n", length(unique(combined_pvals$embedding))))
  cat(sprintf("Embeddings found: %s\n", paste(sort(unique(combined_pvals$embedding)), collapse = ", ")))
  cat(sprintf("Unique CITs: %d\n", length(unique(combined_pvals$cit))))
  cat(sprintf("CITs found: %s\n", paste(sort(unique(combined_pvals$cit)), collapse = ", ")))
  
  # Get available embeddings and CITs
  available_embeddings <- unique(combined_pvals$embedding)
  available_cits <- unique(combined_pvals$cit)
  
} else {
  stop("No p-value data loaded. Check if files exist in the specified directories.")
}

# ============================================================================
# Create Boxplots Per CIT and Condition
# ============================================================================
cat("\n\n=== Creating Boxplots for Each CIT ===\n")

# Define embeddings available for plotting
constant_embeddings_plot <- intersect(constant_embeddings, available_embeddings)
varying_embeddings_plot <- intersect(varying_embeddings, available_embeddings)

cat("Constant embeddings available:", paste(constant_embeddings_plot, collapse = ", "), "\n")
cat("Varying embeddings available:", paste(varying_embeddings_plot, collapse = ", "), "\n")

# Prepare plot data and replace 0 or extremely small p-values
plot_data_pval <- combined_pvals %>%
  mutate(
    p_value = ifelse(p_value < min_pval_plot | p_value == 0, min_pval_plot, p_value),
    n_sample = factor(n_sample, levels = n_samples),
    condition = factor(condition, levels = conditions),
    embedding = factor(embedding),
    cit = factor(cit)
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

# ============================================================================
# Create Individual Boxplots for Each CIT
# ============================================================================

for (cit_name in available_cits) {
  cat(sprintf("\nCreating boxplot for CIT: %s\n", cit_name))
  
  # Filter data for this CIT
  plot_data_single_cit <- plot_data_pval %>%
    filter(cit == cit_name, embedding %in% all_embeddings_plot) %>%
    mutate(embedding = factor(embedding, levels = all_embeddings_plot))
  
  # Check if we have data for this CIT
  if (nrow(plot_data_single_cit) == 0) {
    warning(sprintf("No data available for CIT: %s", cit_name))
    next
  }
  
  # Create boxplot with facets by condition
  p_single_cit <- plot_data_single_cit %>%
    ggplot(aes(x = n_sample, y = p_value, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = "P-value",
      fill = "Embedding",
      title = bquote("CIT P-values for" ~ .(cit_name) ~ "(" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
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
    scale_fill_manual(values = color_palette_pval) +
    scale_y_log10(
      breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.05, 0.1, 0.5, 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = y_limits_pval
    )
  
  # Save plots
  cat("Saving individual CIT boxplots:\n")
  
  # Save as PNG
  png_path_single <- file.path(output_dir, 
                                sprintf("pvalues_embeddings_%s_%d_%d_eps_sigmaY_%s.png", 
                                        cit_name, min(seeds), max(seeds), eps_sigmaY))
  ggsave(
    png_path_single, 
    plot = p_single_cit, 
    width = 16, 
    height = 7, 
    units = "in", 
    dpi = 300
  )
  cat("  Saved PNG:", png_path_single, "\n")
  
  # Save as PDF
  pdf_path_single <- file.path(output_dir, 
                                sprintf("pvalues_embeddings_%s_%d_%d_eps_sigmaY_%s.pdf", 
                                        cit_name, min(seeds), max(seeds), eps_sigmaY))
  ggsave(
    pdf_path_single, 
    plot = p_single_cit, 
    width = 16, 
    height = 7, 
    units = "in"
  )
  cat("  Saved PDF:", pdf_path_single, "\n")
}

cat("\n=== P-value Visualization Complete ===\n")
cat(sprintf("Output directory: %s\n", output_dir))
