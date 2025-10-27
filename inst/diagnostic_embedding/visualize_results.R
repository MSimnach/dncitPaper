# Visualize Diagnostic Results
# This script reads R² results from auto_diagnostic runs across different
# sample sizes, seeds, and eps_sigmaY values, and creates a 4-panel figure.

library(dplyr)
library(ggplot2)
library(data.table)

# ============================================================================
# 1. Data Collection
# ============================================================================
cat("=== Collecting Diagnostic Results ===\n")

# Parameters
n_samples <- c(460, 825, 1100)
seeds <- 1:10
eps_sigmaY_all <- c(0, 0.05, 0.1, 0.5, 1)
base_path <- "/sc/home/marco.simnacher/ukbiobank/data/No_CI"

# Initialize list to store all results
all_results <- list()

# Loop through all combinations and read results
for (n_sample in n_samples) {
  for (seed in seeds) {
    for (eps_sigmaY in eps_sigmaY_all) {
      # Construct path to diagnostic results
      results_path <- file.path(
        base_path, 
        n_sample, 
        seed, 
        paste0("eps_sigmaY=", eps_sigmaY),
        "diagnostic",
        "diagnostic_results.csv"
      )
      
      # Check if file exists
      if (file.exists(results_path)) {
        # Read results
        results_df <- fread(results_path, nThread = 1)
        
        # Add metadata
        results_df$n_sample <- n_sample
        results_df$seed <- seed
        results_df$eps_sigmaY <- eps_sigmaY
        
        # Append to list
        all_results[[length(all_results) + 1]] <- results_df
        
      } else {
        warning("File not found: ", results_path)
      }
    }
  }
}

# Combine all results
if (length(all_results) == 0) {
  stop("No results files found. Please check the paths.")
}

combined_results <- rbindlist(all_results)
cat("Total results loaded:", nrow(combined_results), "\n")
cat("Unique embeddings:", paste(unique(combined_results$embedding_name), collapse = ", "), "\n\n")

# ============================================================================
# 2. Filter and Prepare Data for Plotting
# ============================================================================
cat("=== Preparing Data for Visualization ===\n")

# Filter to only the eps_sigmaY values we want to visualize
plot_data <- combined_results %>%
  filter(eps_sigmaY %in% c(0, 0.1, 0.5, 1)) %>%
  mutate(
    # Convert to factors for proper ordering
    n_sample = factor(n_sample, levels = c(460, 825, 1100)),
    eps_sigmaY = factor(
      eps_sigmaY, 
      levels = c(0, 0.1, 0.5, 1),
      labels = c("eps_sigmaY = 0", "eps_sigmaY = 0.1", 
                 "eps_sigmaY = 0.5", "eps_sigmaY = 1")
    ),
    # Capitalize embedding names for better display
    embedding_name = factor(embedding_name)
  )

cat("Data prepared for plotting:\n")
cat("  - Samples per eps_sigmaY:", table(plot_data$eps_sigmaY), "\n")
cat("  - Samples per n_sample:", table(plot_data$n_sample), "\n")
cat("  - Embeddings:", levels(plot_data$embedding_name), "\n\n")

# ============================================================================
# 3. Create Visualization
# ============================================================================
cat("=== Creating Visualization ===\n")

# Create the plot
p <- ggplot(plot_data, aes(x = n_sample, y = r2_test, fill = embedding_name)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 1) +
  facet_wrap(~ eps_sigmaY, nrow = 2, ncol = 2) +
  labs(
    x = "Sample Size",
    y = expression(R^2 ~ "Test"),
    fill = "Embedding",
    title = "Diagnostic R² Results Across Sample Sizes and Noise Levels"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_fill_brewer(palette = "Set2")

# ============================================================================
# 4. Save Outputs
# ============================================================================
cat("=== Saving Figure ===\n")

# Create output directory
output_dir <- "/sc/home/marco.simnacher/dncitPaper/inst/diagnostic_embedding/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save as PNG
png_path <- file.path(output_dir, "diagnostic_r2_by_sample_size.png")
ggsave(
  png_path, 
  plot = p, 
  width = 12, 
  height = 10, 
  units = "in", 
  dpi = 300
)
cat("Saved PNG:", png_path, "\n")

# Save as PDF
pdf_path <- file.path(output_dir, "diagnostic_r2_by_sample_size.pdf")
ggsave(
  pdf_path, 
  plot = p, 
  width = 12, 
  height = 10, 
  units = "in"
)
cat("Saved PDF:", pdf_path, "\n")

cat("\n=== Visualization Complete ===\n")
cat("Figures saved to:", output_dir, "\n")

# Print summary statistics
cat("\n=== Summary Statistics ===\n")
summary_stats <- plot_data %>%
  group_by(eps_sigmaY, n_sample, embedding_name) %>%
  summarise(
    mean_r2 = mean(r2_test, na.rm = TRUE),
    sd_r2 = sd(r2_test, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
print(summary_stats)

