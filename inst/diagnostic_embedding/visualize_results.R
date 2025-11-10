# Visualize Diagnostic Results
# This script reads R² results from auto_diagnostic runs across different
# sample sizes and seeds for eps_sigmaY = 1, and creates 2x2 grid plots:
# - Left panels: Constant embeddings (fastsurfer, freesurfer, condVAE, medicalnet)
# - Right panels: Trained embeddings (medicalnet_ft, medicalnet_ft_frozen, scratch)
# - Top row: CI condition, Bottom row: No_CI condition
# Creates separate figures for Y diagnostics and age diagnostics

library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(RColorBrewer)

# ============================================================================
# 1. Data Collection
# ============================================================================
cat("=== Collecting Diagnostic Results ===\n")

# Parameters
n_samples <- c(460, 1100, 5000, 10000)#c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)
seeds <- c(301,303, 312:313)#c(51:60, 61:64, 71:74)#1:200
eps_sigmaY_all <- c(1)
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

# Combine all results
if (length(all_results_y) == 0 && length(all_results_age) == 0) {
  stop("No results files found. Please check the paths.")
}

cat("Y Diagnostics:\n")
if (length(all_results_y) > 0) {
  combined_results_y <- rbindlist(all_results_y)
  cat("  Total results loaded:", nrow(combined_results_y), "\n")
  cat("  Unique embeddings:", paste(unique(combined_results_y$embedding_name), collapse = ", "), "\n")
} else {
  warning("No Y diagnostic results found.")
  combined_results_y <- data.table()
}

cat("\nAge Diagnostics:\n")
if (length(all_results_age) > 0) {
  combined_results_age <- rbindlist(all_results_age)
  cat("  Total results loaded:", nrow(combined_results_age), "\n")
  cat("  Unique embeddings:", paste(unique(combined_results_age$embedding_name), collapse = ", "), "\n")
} else {
  warning("No age diagnostic results found.")
  combined_results_age <- data.table()
}
cat("\n")

# ============================================================================
# 2. Filter and Prepare Data for Plotting
# ============================================================================
cat("=== Preparing Data for Visualization ===\n")

# Define embedding types
constant_embeddings <- c("fastsurfer", "freesurfer", "condVAE", "medicalnet")
varying_embeddings <- c("medicalnet_ft", "medicalnet_ft_frozen", "scratch")

# Function to prepare plot data
prepare_plot_data <- function(combined_results) {
  if (nrow(combined_results) == 0) {
    return(NULL)
  }
  
  plot_data_all <- combined_results %>%
    filter(eps_sigmaY == 1) %>%
    mutate(
      # Convert to factors for proper ordering
      n_sample = factor(n_sample, levels = c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
      condition = factor(condition, levels = c("CI", "No_CI")),
      # Capitalize embedding names for better display
      embedding_name = factor(embedding_name)
    )
  
  # Split into constant and varying embeddings
  # For constant embeddings, use only the largest sample size
  plot_data_constant <- plot_data_all %>%
    filter(embedding_name %in% constant_embeddings, n_sample == "10000") %>%
    mutate(embedding_name = factor(embedding_name, levels = constant_embeddings))
  
  plot_data_varying <- plot_data_all %>%
    filter(embedding_name %in% varying_embeddings) %>%
    mutate(embedding_name = factor(embedding_name, levels = varying_embeddings))
  
  return(list(
    all = plot_data_all,
    constant = plot_data_constant,
    varying = plot_data_varying
  ))
}

# Prepare data for Y diagnostics
cat("\nPreparing Y diagnostic data:\n")
plot_data_y <- prepare_plot_data(combined_results_y)
if (!is.null(plot_data_y)) {
  cat("  - Y constant embeddings:", nrow(plot_data_y$constant), "rows\n")
  cat("  - Y varying embeddings:", nrow(plot_data_y$varying), "rows\n")
}

# Prepare data for age diagnostics
cat("\nPreparing age diagnostic data:\n")
plot_data_age <- prepare_plot_data(combined_results_age)
if (!is.null(plot_data_age)) {
  cat("  - Age constant embeddings:", nrow(plot_data_age$constant), "rows\n")
  cat("  - Age varying embeddings:", nrow(plot_data_age$varying), "rows\n")
}

cat("\nData summary:\n")
cat("  - eps_sigmaY: 1\n")
cat("  - Conditions: CI, No_CI\n")
cat("  - Constant embeddings (n=10000):", paste(constant_embeddings, collapse=", "), "\n")
cat("  - Varying embeddings:", paste(varying_embeddings, collapse=", "), "\n\n")

# ============================================================================
# 3. Create Visualization
# ============================================================================
cat("=== Creating Visualization ===\n")

# Create color palette for all embeddings
all_embeddings <- c(constant_embeddings, varying_embeddings)
color_palette <- setNames(
  RColorBrewer::brewer.pal(length(all_embeddings), "Set2"),
  all_embeddings
)

# Function to create 2x2 grid plot for a diagnostic type
create_diagnostic_plot <- function(plot_data, y_limits, diagnostic_name) {
  if (is.null(plot_data)) {
    warning(paste("No data available for", diagnostic_name, "diagnostics"))
    return(NULL)
  }
  
  # Panel 1: Top Left - CI, Constant embeddings
  p1 <- plot_data$constant %>%
    filter(condition == "CI") %>%
    ggplot(aes(x = embedding_name, y = r2_test, fill = embedding_name)) +
    geom_boxplot(outlier.size = 1) +
    labs(
      x = "Embedding",
      y = expression(R^2 ~ "Test"),
      title = "CI - Constant Embeddings"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = color_palette) +
    coord_cartesian(ylim = y_limits)
  
  # Panel 2: Bottom Left - No_CI, Constant embeddings
  p2 <- plot_data$constant %>%
    filter(condition == "No_CI") %>%
    ggplot(aes(x = embedding_name, y = r2_test, fill = embedding_name)) +
    geom_boxplot(outlier.size = 1) +
    labs(
      x = "Embedding",
      y = expression(R^2 ~ "Test"),
      title = "No_CI - Constant Embeddings"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = color_palette) +
    coord_cartesian(ylim = y_limits)
  
  # Panel 3: Top Right - CI, Varying embeddings across sample sizes
  p3 <- plot_data$varying %>%
    filter(condition == "CI") %>%
    ggplot(aes(x = n_sample, y = r2_test, fill = embedding_name)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 1) +
    labs(
      x = "Sample Size",
      y = expression(R^2 ~ "Test"),
      title = "CI - Trained Embeddings"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = color_palette) +
    coord_cartesian(ylim = y_limits)
  
  # Panel 4: Bottom Right - No_CI, Varying embeddings across sample sizes
  p4 <- plot_data$varying %>%
    filter(condition == "No_CI") %>%
    ggplot(aes(x = n_sample, y = r2_test, fill = embedding_name)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 1) +
    labs(
      x = "Sample Size",
      y = expression(R^2 ~ "Test"),
      title = "No_CI - Trained Embeddings"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = color_palette) +
    coord_cartesian(ylim = y_limits)
  
  # Combine plots in 2x2 grid with shared legend at bottom
  # Set width ratio of 1:3 for left:right panels
  top_row <- (p1 | p3) + plot_layout(widths = c(1, 3))
  bottom_row <- (p2 | p4) + plot_layout(widths = c(1, 3))
  p_combined <- top_row / bottom_row + 
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.title = element_text(face = "bold"))
  
  # Add overall title
  title_text <- if (diagnostic_name == "Y") {
    expression("Diagnostic R² Results: CI vs No_CI (" * epsilon[sigma[Y]] * " = 1)")
  } else {
    expression("Age Diagnostic R² Results: CI vs No_CI (" * epsilon[sigma[Y]] * " = 1)")
  }
  
  p <- p_combined + 
    plot_annotation(
      title = title_text,
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
  
  return(p)
}

# Calculate y-axis limits for each diagnostic type
cat("\nCalculating y-axis limits:\n")
if (!is.null(plot_data_y)) {
  y_limits_y <- c(
    min(plot_data_y$all$r2_test, na.rm = TRUE),
    max(plot_data_y$all$r2_test, na.rm = TRUE)
  )
  cat("  Y diagnostics: [", y_limits_y[1], ", ", y_limits_y[2], "]\n")
} else {
  y_limits_y <- NULL
}

if (!is.null(plot_data_age)) {
  y_limits_age <- c(
    min(plot_data_age$all$r2_test, na.rm = TRUE),
    max(plot_data_age$all$r2_test, na.rm = TRUE)
  )
  cat("  Age diagnostics: [", y_limits_age[1], ", ", y_limits_age[2], "]\n")
} else {
  y_limits_age <- NULL
}
cat("\n")

# Create plots for Y diagnostics
cat("Creating Y diagnostic plot...\n")
p_y <- create_diagnostic_plot(plot_data_y, y_limits_y, "Y")

# Create plots for age diagnostics
cat("Creating age diagnostic plot...\n")
p_age <- create_diagnostic_plot(plot_data_age, y_limits_age, "Age")

# ============================================================================
# 4. Save Outputs
# ============================================================================
cat("\n=== Saving Figures ===\n")

# Create output directory
output_dir <- "/sc/home/marco.simnacher/dncitPaper/inst/diagnostic_embedding/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save Y diagnostics plot
if (!is.null(p_y)) {
  cat("\nSaving Y diagnostic plots:\n")
  # Save as PNG
  png_path_y <- file.path(output_dir, "diagnostic_r2_by_sample_size.png")
  ggsave(
    png_path_y, 
    plot = p_y, 
    width = 16, 
    height = 10, 
    units = "in", 
    dpi = 300
  )
  cat("  Saved PNG:", png_path_y, "\n")
  
  # Save as PDF
  pdf_path_y <- file.path(output_dir, "diagnostic_r2_by_sample_size.pdf")
  ggsave(
    pdf_path_y, 
    plot = p_y, 
    width = 16, 
    height = 10, 
    units = "in"
  )
  cat("  Saved PDF:", pdf_path_y, "\n")
}

# Save age diagnostics plot
if (!is.null(p_age)) {
  cat("\nSaving age diagnostic plots:\n")
  # Save as PNG
  png_path_age <- file.path(output_dir, "diagnostic_r2_by_sample_size_age.png")
  ggsave(
    png_path_age, 
    plot = p_age, 
    width = 16, 
    height = 10, 
    units = "in", 
    dpi = 300
  )
  cat("  Saved PNG:", png_path_age, "\n")
  
  # Save as PDF
  pdf_path_age <- file.path(output_dir, "diagnostic_r2_by_sample_size_age.pdf")
  ggsave(
    pdf_path_age, 
    plot = p_age, 
    width = 16, 
    height = 10, 
    units = "in"
  )
  cat("  Saved PDF:", pdf_path_age, "\n")
}

cat("\n=== Visualization Complete ===\n")
cat("Figures saved to:", output_dir, "\n")

# Print summary statistics
cat("\n=== Summary Statistics (eps_sigmaY = 1) ===\n")

if (!is.null(plot_data_y)) {
  cat("\n--- Y Diagnostics ---\n")
  
  cat("\nConstant Embeddings (n=10000):\n")
  summary_stats_y_constant <- plot_data_y$constant %>%
    group_by(condition, embedding_name) %>%
    summarise(
      mean_r2 = mean(r2_test, na.rm = TRUE),
      sd_r2 = sd(r2_test, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  print(summary_stats_y_constant)
  
  cat("\nVarying Embeddings (across sample sizes):\n")
  summary_stats_y_varying <- plot_data_y$varying %>%
    group_by(condition, n_sample, embedding_name) %>%
    summarise(
      mean_r2 = mean(r2_test, na.rm = TRUE),
      sd_r2 = sd(r2_test, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  print(summary_stats_y_varying)
}

if (!is.null(plot_data_age)) {
  cat("\n--- Age Diagnostics ---\n")
  
  cat("\nConstant Embeddings (n=10000):\n")
  summary_stats_age_constant <- plot_data_age$constant %>%
    group_by(condition, embedding_name) %>%
    summarise(
      mean_r2 = mean(r2_test, na.rm = TRUE),
      sd_r2 = sd(r2_test, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  print(summary_stats_age_constant)
  
  cat("\nVarying Embeddings (across sample sizes):\n")
  summary_stats_age_varying <- plot_data_age$varying %>%
    group_by(condition, n_sample, embedding_name) %>%
    summarise(
      mean_r2 = mean(r2_test, na.rm = TRUE),
      sd_r2 = sd(r2_test, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  print(summary_stats_age_varying)
}

