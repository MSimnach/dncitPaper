### P-value Extraction and Boxplot Visualization
# Extract CIT p-values from Results/ directory and create boxplots
# showing p-values by CIT, embedding, condition, and sample size

# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(paletteer)

# ============================================================================
# Configuration Variables
# ============================================================================

# Seeds to analyze
seeds <- c(1:100)

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

# Function to load rejection rates from CSV file
load_rejection_rates <- function(csv_path) {
  if (!file.exists(csv_path)) {
    return(NULL)
  }
  rr_df <- fread(csv_path, nThread = 1)
  colnames(rr_df) <- c("n_sample", "rejection_rate")
  return(rr_df)
}

# ============================================================================
# Load P-values from CSV Files
# ============================================================================

cat("\n=== Loading P-values from CIT Tests ===\n")

# Construct seed directory name
seed_dir_name <- sprintf("seeds_%d_%d", min(seeds), max(seeds))

# Initialize list to store p-value data
pval_data_list <- list()

# Initialize list to store rejection rate data
rr_data_list <- list()

for (cond in conditions) {
  pval_dir <- file.path(results_base_path, cond, "p-values", seed_dir_name)
  rr_dir <- file.path(results_base_path, cond, "rejection_rates", seed_dir_name)
  
  # Check if directory exists
  if (!dir.exists(pval_dir)) {
    warning(sprintf("Directory not found: %s", pval_dir))
    next
  }
  
  # Get all CSV files in the directory that match the specified eps_sigmaY
  # Pattern: 1_0_[eps_sigmaY]_0_...csv
  pattern_with_eps <- sprintf("^1_0_%g_0_.*\\.csv$", eps_sigmaY)
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
    
    # Load corresponding rejection rates
    rr_file <- file.path(rr_dir, filename)
    if (file.exists(rr_file)) {
      rr_df <- load_rejection_rates(rr_file)
      if (!is.null(rr_df)) {
        # Add metadata
        rr_df$embedding <- embedding_name
        rr_df$condition <- cond
        rr_df$cit <- cit_name
        
        # Store in list
        rr_data_list[[length(rr_data_list) + 1]] <- rr_df
      }
    } else {
      warning(sprintf("Rejection rate file not found: %s", rr_file))
    }
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

# Combine all rejection rate data
if (length(rr_data_list) > 0) {
  combined_rr <- rbindlist(rr_data_list)
  
  cat(sprintf("\n=== Rejection Rate Summary ===\n"))
  cat(sprintf("Total rejection rate records loaded: %d\n", nrow(combined_rr)))
  cat(sprintf("Unique embeddings in RR data: %d\n", length(unique(combined_rr$embedding))))
  cat(sprintf("Unique CITs in RR data: %d\n", length(unique(combined_rr$cit))))
  
} else {
  warning("No rejection rate data loaded. Skipping rejection rate plots.")
  combined_rr <- NULL
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

# Conditionally include medicalnet_ft based on eps_sigmaY
if (eps_sigmaY == 0.5 && "medicalnet_ft" %in% available_embeddings) {
  embeddings_to_plot <- c(constant_embeddings_plot, "medicalnet_ft")
  if ("scratch" %in% varying_embeddings_plot) {
    embeddings_to_plot <- c(embeddings_to_plot, "scratch")
  }
  cat("Including medicalnet_ft for eps_sigmaY = 0.5\n")
} else {
  embeddings_to_plot <- c(constant_embeddings_plot, 
                          setdiff(varying_embeddings_plot, "medicalnet_ft"))
  cat("Excluding medicalnet_ft (eps_sigmaY != 0.5)\n")
}

# Define embedding order: fastsurfer, freesurfer, condVAE, medicalnet, medicalnet_ft, scratch
embedding_order <- c("fastsurfer", "freesurfer", "condVAE", "medicalnet", "scratch", "medicalnet_ft")
all_embeddings_plot <- intersect(embedding_order, embeddings_to_plot)

cat("Embeddings to plot:", paste(all_embeddings_plot, collapse = ", "), "\n")

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

# Create color palette - maintain consistent colors
# Colors 1-4: fastsurfer, freesurfer, condVAE, medicalnet (always same)
# Color 5: medicalnet_ft (when included)
# Color 6: scratch (when included)
palet_discrete <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
color_palette_pval <- setNames(
  palet_discrete[1:length(all_embeddings_plot)],
  all_embeddings_display
)

# Prepare plot data and replace 0 or extremely small p-values
plot_data_pval <- combined_pvals %>%
  mutate(
    p_value = ifelse(p_value < min_pval_plot | p_value == 0, min_pval_plot, p_value),
    n_sample = factor(n_sample, levels = n_samples),
    condition = factor(condition, levels = conditions),
    embedding = factor(embedding, levels = all_embeddings_plot, labels = all_embeddings_display),
    cit = factor(cit)
  ) %>%
  mutate(
    condition = factor(condition, levels = conditions, labels = c("T1E", "Power"))
  )

# Report how many p-values were adjusted
n_adjusted <- sum(combined_pvals$p_value < min_pval_plot | combined_pvals$p_value == 0, na.rm = TRUE)
if (n_adjusted > 0) {
  cat(sprintf("Note: %d p-values (%.1f%%) were < %g or zero, set to %g for plotting\n", 
              n_adjusted, 100*n_adjusted/nrow(combined_pvals), min_pval_plot, min_pval_plot))
}

# Calculate y-axis limits (log scale) - use adjusted data so no Inf values
y_limits_pval <- c(min_pval_plot, 1)

# ============================================================================
# Create Individual Boxplots for Each CIT
# ============================================================================

for (cit_name in available_cits) {
  cat(sprintf("\nCreating boxplot for CIT: %s\n", cit_name))
  
  # Filter data for this CIT
  plot_data_single_cit <- plot_data_pval %>%
    filter(cit == cit_name, embedding %in% all_embeddings_display) %>%
    mutate(embedding = factor(embedding, levels = all_embeddings_display))
  
  # Check if we have data for this CIT
  if (nrow(plot_data_single_cit) == 0) {
    warning(sprintf("No data available for CIT: %s", cit_name))
    next
  }
  
  # Create boxplot with facets by condition
  p_single_cit <- plot_data_single_cit %>%
    ggplot(aes(x = n_sample, y = p_value, fill = embedding)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
    geom_hline(yintercept = 0.05, color = "black", alpha = 0.6) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = "P-value",
      fill = "Embedding"#,
      # title = bquote("CIT P-values for" ~ .(cit_name) ~ "(" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      legend.position = "right",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(size = 20)
    ) +
    scale_fill_manual(values = color_palette_pval) +
    scale_x_discrete(breaks = c("256", "460", "825", "1100", "5000")) +
    scale_y_log10(
      breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 1),
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

# ============================================================================
# Create Combined Rejection Rate Line Plot
# ============================================================================

if (!is.null(combined_rr)) {
  cat("\n\n=== Creating Combined Rejection Rate Line Plot ===\n")
  
  # Prepare rejection rate plot data
  plot_data_rr <- combined_rr %>%
    filter(cit %in% available_cits, embedding %in% all_embeddings_plot) %>%
    mutate(
      n_sample = as.numeric(as.character(n_sample)),
      embedding = factor(embedding, levels = all_embeddings_plot, labels = all_embeddings_display),
      cit = factor(cit, levels = available_cits),
      condition = factor(condition, levels = conditions)
    ) %>%
    mutate(
      condition = factor(condition, levels = conditions, labels = c("T1E", "Power"))
    ) %>%
    mutate(
      # Create equidistant position mapping
      n_sample_pos = match(n_sample, sort(unique(n_sample)))
    )
  
  # Create color palette for embeddings (same as plot_results_nested_loop.R)
  palet_discrete <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
  embedding_colors <- setNames(
    palet_discrete[1:length(all_embeddings_plot)],
    all_embeddings_display
  )
  
  # Create shape and linetype mappings for CITs
  cit_shapes <- setNames(
    c(16, 17, 15, 18, 8, 4)[1:length(available_cits)], 
    available_cits
  )
  cit_linetypes <- setNames(
    c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")[1:length(available_cits)],
    available_cits
  )
  
  # Create line plot with equidistant x-axis
  # Get unique sample sizes in order for labels
  sample_size_labels <- sort(unique(plot_data_rr$n_sample))
  sample_size_positions <- 1:length(sample_size_labels)
  
  p_rejection_rates_combined <- plot_data_rr %>%
    ggplot(aes(x = n_sample_pos, y = rejection_rate, 
               color = embedding, shape = cit, linetype = cit,
               group = interaction(cit, embedding))) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    #geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
    facet_wrap(~ condition, ncol = 2) +
    labs(
      x = "Sample Size",
      y = "Rejection Rate",
      color = "Embedding",
      shape = "CIT",
      linetype = "CIT"#,
      # title = bquote("Rejection Rates Across CITs and Embeddings (" * epsilon[sigma[Y]] * " = " * .(eps_sigmaY) * ")")
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "lightgray"),
      strip.text = element_text(size = 20)
    ) +
    scale_color_manual(values = embedding_colors) +
    scale_shape_manual(values = cit_shapes) +
    scale_linetype_manual(values = cit_linetypes) +
    scale_x_continuous(breaks = sample_size_positions, labels = sample_size_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))
  
  # Save combined rejection rate plot
  cat("\nSaving combined rejection rate plot:\n")
  
  # Save as PNG
  png_path_combined_rr <- file.path(output_dir, 
                                     sprintf("rejection_rates_combined_%d_%d_eps_sigmaY_%s.png", 
                                             min(seeds), max(seeds), eps_sigmaY))
  ggsave(
    png_path_combined_rr, 
    plot = p_rejection_rates_combined, 
    width = 18, 
    height = 10, 
    units = "in", 
    dpi = 300
  )
  cat("  Saved PNG:", png_path_combined_rr, "\n")
  
  # Save as PDF
  pdf_path_combined_rr <- file.path(output_dir, 
                                     sprintf("rejection_rates_combined_%d_%d_eps_sigmaY_%s.pdf", 
                                             min(seeds), max(seeds), eps_sigmaY))
  ggsave(
    pdf_path_combined_rr, 
    plot = p_rejection_rates_combined, 
    width = 18, 
    height = 10, 
    units = "in"
  )
  cat("  Saved PDF:", pdf_path_combined_rr, "\n")
  
  cat("\n=== Rejection Rate Visualization Complete ===\n")
}

cat("\n=== P-value Visualization Complete ===\n")
cat(sprintf("Output directory: %s\n", output_dir))

# ============================================================================
# Create Combined 3-Panel Plots (Rejection Rates + PCM + RCOT p-values)
# ============================================================================

cat("\n\n=== Creating Combined 3-Panel Plots ===\n")

# Filter data for specific CITs
cit_pcm <- "comets_pcm"
cit_rcot <- "RCOT_1"

# Check if both CITs exist in the data
if (cit_pcm %in% available_cits && cit_rcot %in% available_cits && !is.null(combined_rr)) {
  
  # Filter p-value data for PCM and RCOT
  pval_data_pcm <- plot_data_pval %>%
    filter(cit == cit_pcm, embedding %in% all_embeddings_display) %>%
    mutate(embedding = factor(embedding, levels = all_embeddings_display))
  
  pval_data_rcot <- plot_data_pval %>%
    filter(cit == cit_rcot, embedding %in% all_embeddings_display) %>%
    mutate(embedding = factor(embedding, levels = all_embeddings_display))
  
  # Filter rejection rate data (all CITs)
  rr_data_combined <- plot_data_rr
  
  # Rename CIT labels
  pval_data_pcm <- pval_data_pcm %>%
    mutate(cit = factor(cit, levels = c(cit_pcm), labels = c("PCM")))
  
  pval_data_rcot <- pval_data_rcot %>%
    mutate(cit = factor(cit, levels = c(cit_rcot), labels = c("RCoT")))
  
  rr_data_combined <- rr_data_combined %>%
    mutate(cit = factor(cit, 
                        levels = available_cits,
                        labels = ifelse(available_cits == cit_pcm, "PCM",
                                      ifelse(available_cits == cit_rcot, "RCoT", 
                                            as.character(available_cits)))))
  
  # Create separate combined plots for T1E and Power
  for (cond in c("T1E", "Power")) {
    cat(sprintf("\nCreating 3-panel plot for condition: %s\n", cond))
    
    # Filter data for this condition
    pval_pcm_cond <- pval_data_pcm %>% filter(condition == cond)
    pval_rcot_cond <- pval_data_rcot %>% filter(condition == cond)
    rr_cond <- rr_data_combined %>% filter(condition == cond)
    
    # Get unique sample sizes and positions from the actual data
    pos_label_map <- rr_cond %>%
      select(n_sample_pos, n_sample) %>%
      distinct() %>%
      arrange(n_sample_pos)
    sample_size_positions <- pos_label_map$n_sample_pos
    sample_size_labels <- pos_label_map$n_sample
    
    # Update CIT colors/shapes for renamed labels
    cit_labels_renamed <- unique(rr_cond$cit)
    cit_shapes_renamed <- setNames(
      c(16, 17, 15, 18, 8, 4)[1:length(cit_labels_renamed)], 
      cit_labels_renamed
    )
    cit_linetypes_renamed <- setNames(
      c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")[1:length(cit_labels_renamed)],
      cit_labels_renamed
    )
    
    # Panel 1: Rejection Rates (top)
    p_rr_panel <- rr_cond %>%
      ggplot(aes(x = n_sample_pos, y = rejection_rate, 
                 color = embedding, shape = cit, linetype = cit,
                 group = interaction(cit, embedding))) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      labs(
        y = "Rejection Rate",
        color = "Embedding",
        shape = "CIT",
        linetype = "CIT"
      ) +
      theme_bw(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
      ) +
      scale_color_manual(values = embedding_colors) +
      scale_shape_manual(values = cit_shapes_renamed) +
      scale_linetype_manual(values = cit_linetypes_renamed) +
      scale_x_continuous(breaks = sample_size_positions, labels = sample_size_labels)
    
    # Add y-axis scale based on condition
    if (cond == "T1E") {
      p_rr_panel <- p_rr_panel +
        scale_y_continuous(limits = c(0, max(rr_cond$rejection_rate, na.rm = TRUE) * 1.05), 
                          breaks = scales::pretty_breaks(n = 10)) +
                          geom_hline(yintercept = 0.05, color = "black", alpha = 0.6) 
    } else {
      p_rr_panel <- p_rr_panel +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))
    }
    
    # Panel 2: PCM P-values (middle)
    p_pcm_panel <- pval_pcm_cond %>%
      ggplot(aes(x = n_sample, y = p_value, fill = embedding)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
      geom_hline(yintercept = 0.05, color = "black", alpha = 0.6) +
      labs(
        y = "P-value (PCM)",
        fill = "Embedding"
      ) +
      theme_bw(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
      ) +
      scale_fill_manual(values = color_palette_pval) +
      scale_x_discrete(breaks = c("256", "460", "825", "1100", "5000"))
    
    # Add y-axis scale based on condition for PCM panel
    if (cond == "T1E") {
      p_pcm_panel <- p_pcm_panel +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))
    } else {
      p_pcm_panel <- p_pcm_panel +
        scale_y_log10(
          breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 1),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          limits = y_limits_pval
        )
    }
    p_pcm_panel <- p_pcm_panel + guides(fill = "none")
    
    # Panel 3: RCOT P-values (bottom)
    p_rcot_panel <- pval_rcot_cond %>%
      ggplot(aes(x = n_sample, y = p_value, fill = embedding)) +
      geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
      geom_hline(yintercept = 0.05, color = "black", alpha = 0.6) +
      labs(
        x = "Sample Size",
        y = "P-value (RCoT)",
        fill = "Embedding"
      ) +
      theme_bw(base_size = 18) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16)
      ) +
      scale_fill_manual(values = color_palette_pval) +
      scale_x_discrete(breaks = c("256", "460", "825", "1100", "5000"))
    
    # Add y-axis scale based on condition for RCoT panel
    if (cond == "T1E") {
      p_rcot_panel <- p_rcot_panel +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))
    } else {
      p_rcot_panel <- p_rcot_panel +
        scale_y_log10(
          breaks = c(1e-16, 1e-12, 1e-9, 1e-6, 1e-3, 0.01, 0.1, 1),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          limits = y_limits_pval
        )
    }
    p_rcot_panel <- p_rcot_panel + guides(fill = "none")
    
    # Combine plots vertically with shared legend
    p_combined_3panel <- p_rr_panel / p_pcm_panel / p_rcot_panel +
      plot_layout(ncol = 1, heights = c(1, 1, 1), guides = "collect") &
      theme(legend.position = "right")
    
    # Save combined 3-panel plot
    cat(sprintf("Saving 3-panel plot for %s:\n", cond))
    
    # Save as PNG
    png_path_3panel <- file.path(output_dir, 
                                  sprintf("combined_3panel_%s_%d_%d_eps_sigmaY_%s.png", 
                                          cond, min(seeds), max(seeds), eps_sigmaY))
    ggsave(
      png_path_3panel, 
      plot = p_combined_3panel, 
      width = 16, 
      height = 18, 
      units = "in", 
      dpi = 300
    )
    cat("  Saved PNG:", png_path_3panel, "\n")
    
    # Save as PDF
    pdf_path_3panel <- file.path(output_dir, 
                                  sprintf("combined_3panel_%s_%d_%d_eps_sigmaY_%s.pdf", 
                                          cond, min(seeds), max(seeds), eps_sigmaY))
    ggsave(
      pdf_path_3panel, 
      plot = p_combined_3panel, 
      width = 16, 
      height = 18, 
      units = "in"
    )
    cat("  Saved PDF:", pdf_path_3panel, "\n")
    
  }
  
  cat("\n=== Combined 3-Panel Plots Complete ===\n")
  
} else {
  warning("Cannot create combined 3-panel plots: PCM and/or RCOT not found in data, or rejection rate data missing.")
}

# ============================================================================
# Create QQ Plots for eps_sigmaY = 0.5
# ============================================================================

if (eps_sigmaY == 0.5) {
  cat("\n\n=== Creating QQ Plots ===\n")
  
  # Load required library for grid arrangement
  library(gridExtra)
  library(grid)
  
  # CITs to analyze
  cit_pcm <- "comets_pcm"
  cit_rcot <- "RCOT_1"
  
  # Check if both CITs exist
  if (cit_pcm %in% available_cits && cit_rcot %in% available_cits) {
    
    # Prepare data for QQ plots - use plot_data_pval which has renamed conditions
    qq_data <- plot_data_pval %>%
      filter(embedding %in% all_embeddings_display,
             cit %in% c(cit_pcm, cit_rcot))
    
    # Create QQ plots for each CIT and each condition
    for (cit_qq in c(cit_pcm, cit_rcot)) {
      cit_label <- ifelse(cit_qq == cit_pcm, "PCM", "RCoT")
      
      for (cond in c("T1E", "Power")) {
        cat(sprintf("\nCreating QQ plot for %s - %s\n", cit_label, cond))
        
        # Filter data for this CIT and condition
        qq_data_cit <- qq_data %>% 
          filter(cit == cit_qq, condition == cond)
        
        # Number of embeddings and sample sizes
        n_embeddings <- length(all_embeddings_display)
        n_rows <- length(n_samples)
        
        # List to store all plots (organized by rows first, then columns)
        # Layout: Columns = embeddings, Rows = sample sizes (matching original)
        plot_list <- list()
        
        # Loop over each sample size (rows)
        for (row_idx in seq_along(n_samples)) {
          n_s <- n_samples[row_idx]
          
          # Loop over each embedding (columns)
          for (col_idx in seq_along(all_embeddings_display)) {
            emb_display <- all_embeddings_display[col_idx]
            emb_internal <- all_embeddings_plot[col_idx]
            
            # Extract p-values for this embedding and sample size (use display name)
            p_vals <- qq_data_cit %>%
              filter(embedding == emb_display, n_sample == n_s) %>%
              pull(p_value)
            
            # Remove NA values
            p_vals <- p_vals[!is.na(p_vals)]
            
            if (length(p_vals) > 0) {
              # Calculate KS test statistic
              ks_stat <- round(ks.test(p_vals, "punif")$statistic, 3)
              
              # Create data frame for plotting
              plot_df <- data.frame(
                Observed = sort(p_vals),
                Theoretical = qunif(ppoints(length(p_vals)))
              )
              
              # Create the QQ plot
              p <- ggplot(plot_df, aes(sample = Observed)) +
                stat_qq(distribution = stats::qunif, color = color_palette_pval[emb_display]) +
                geom_abline(slope = 1, intercept = 0, color = "black") +
                xlab(ks_stat) +
                ylab(n_s) +
                theme_minimal() +
                coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
              
              # Base theme
              p <- p + theme(
                plot.title = element_blank(),
                text = element_text(size = 15),
                axis.title.y = element_blank()
              )
              
              # Visibility logic matching original implementation
              # First column: show y-axis with sample size label
              if (col_idx == 1) {
                p <- p + theme(
                  axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 5, b = 0, l = 0))
                )
              } else {
                # Other columns: remove y-axis
                p <- p + theme(
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank()
                )
              }
              
              # Last row: show x-axis with KS statistic
              if (row_idx != n_rows) {
                # Not last row: remove x-axis text and ticks
                p <- p + theme(
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()
                )
              }
              
              plot_list[[length(plot_list) + 1]] <- p
            } else {
              # Create empty plot if no data
              p <- ggplot() +
                theme_void() +
                annotate("text", x = 0.5, y = 0.5, label = "No data", size = 5, hjust = 0.5)
              plot_list[[length(plot_list) + 1]] <- p
            }
          }
        }
        
        # Create column headers (embedding display names)
        col_headers <- lapply(all_embeddings_display, function(emb) {
          textGrob(emb, gp = gpar(fontsize = 15), just = "center")
        })
        
        # Combine headers and plots
        grob_list <- c(col_headers, plot_list)
        
        # Create layout matrix matching original:
        # First row: column headers (embeddings)
        # Subsequent rows: plots (sample sizes x embeddings)
        layout_matrix <- rbind(
          matrix(1:n_embeddings, nrow = 1, ncol = n_embeddings),
          matrix((n_embeddings + 1):(n_embeddings + n_embeddings * n_rows), 
                 nrow = n_rows, ncol = n_embeddings, byrow = TRUE)
        )
        
        # Arrange all plots into a grid
        grid_qq <- gridExtra::grid.arrange(
          grobs = grob_list,
          layout_matrix = layout_matrix,
          bottom = textGrob("Theoretical Quantiles", gp = gpar(fontsize = 20)),
          left = textGrob("Sample Quantiles", gp = gpar(fontsize = 20), rot = 90),
          heights = unit.c(unit(1, "lines"), unit(1.1, "null"), rep(unit(1, "null"), n_rows - 1))
        )
        
        # Save QQ plot
        cat(sprintf("Saving QQ plot for %s - %s:\n", cit_label, cond))
        
        # Save as PNG
        png_path_qq <- file.path(output_dir, 
                                  sprintf("qq_plot_%s_%s_%d_%d_eps_sigmaY_%s.png", 
                                          cit_label, cond, min(seeds), max(seeds), eps_sigmaY))
        ggsave(
          png_path_qq, 
          plot = grid_qq, 
          width = n_embeddings * 2.5, 
          height = n_rows * 2.5, 
          units = "in", 
          dpi = 300
        )
        cat("  Saved PNG:", png_path_qq, "\n")
        
        # Save as PDF
        pdf_path_qq <- file.path(output_dir, 
                                  sprintf("qq_plot_%s_%s_%d_%d_eps_sigmaY_%s.pdf", 
                                          cit_label, cond, min(seeds), max(seeds), eps_sigmaY))
        ggsave(
          pdf_path_qq, 
          plot = grid_qq, 
          width = n_embeddings * 2.5, 
          height = n_rows * 2.5, 
          units = "in"
        )
        cat("  Saved PDF:", pdf_path_qq, "\n")
      }
    }
    
    cat("\n=== QQ Plots Complete ===\n")
    
  } else {
    warning("Cannot create QQ plots: PCM and/or RCOT not found in data.")
  }
} else {
  cat("\n=== Skipping QQ Plots (only created for eps_sigmaY = 0.5) ===\n")
}

# ============================================================================
# Create Runtime Boxplots for Scratch vs MedicalNet-ft
# ============================================================================

if (eps_sigmaY == 0.5) {
  cat("\n\n=== Creating Runtime Boxplots ===\n")
  
  # Define runtime embeddings to compare
  runtime_embeddings <- c("scratch", "medicalnet_ft")
  
  # CIT method for runtime data
  cit_runtime <- "comets_pcm"
  
  # Initialize list to store runtime data
  runtime_data_list <- list()
  
  # Load runtime data from both CI and No_CI directories
  for (cond in conditions) {
    # Construct path to runtime directory
    runtime_dir <- file.path(results_base_path, cond, "runtime_embedding", seed_dir_name)
    
    # Check if directory exists
    if (!dir.exists(runtime_dir)) {
      warning(sprintf("Runtime directory not found: %s", runtime_dir))
      next
    }
    
    cat(sprintf("\nLoading runtime data from: %s\n", runtime_dir))
    
    # Load runtime data for each embedding
    for (emb in runtime_embeddings) {
      # Construct filename
      runtime_file <- sprintf("1_0_%g_0_fastsurfer_%s_ukb_z1_squared_%s.csv", 
                               eps_sigmaY, emb, cit_runtime)
      runtime_path <- file.path(runtime_dir, runtime_file)
      
      if (file.exists(runtime_path)) {
        cat(sprintf("  Loading: %s (condition: %s)\n", runtime_file, cond))
        
        # Read the CSV
        runtime_df <- fread(runtime_path, nThread = 1)
        
        # Drop the first column (row indices)
        runtime_df <- runtime_df[, -1]
        
        # Get column names for sample sizes (V1, V2, V3, V4, V5)
        sample_size_cols <- names(runtime_df)[1:min(length(n_samples), ncol(runtime_df))]
        
        # Add seed information based on row number
        runtime_df$seed_idx <- 1:nrow(runtime_df)
        runtime_df$seed <- seeds[runtime_df$seed_idx]
        
        # Reshape to long format
        runtime_long <- melt(runtime_df, 
                            id.vars = c("seed", "seed_idx"),
                            measure.vars = sample_size_cols,
                            variable.name = "sample_size_col",
                            value.name = "runtime")
        
        # Map V1-V5 to actual sample sizes
        runtime_long$n_sample <- n_samples[as.integer(gsub("V", "", runtime_long$sample_size_col))]
        
        # Add embedding and condition metadata
        runtime_long$embedding <- emb
        runtime_long$condition <- cond
        
        # Store in list
        runtime_data_list[[length(runtime_data_list) + 1]] <- runtime_long[, .(seed, n_sample, embedding, condition, runtime)]
        
      } else {
        warning(sprintf("Runtime file not found: %s", runtime_path))
      }
    }
  }
    
    # Combine all runtime data
    if (length(runtime_data_list) > 0) {
      combined_runtime <- rbindlist(runtime_data_list)
      
      cat(sprintf("\n=== Runtime Data Summary ===\n"))
      cat(sprintf("Total runtime records loaded: %d\n", nrow(combined_runtime)))
      cat(sprintf("Embeddings: %s\n", paste(unique(combined_runtime$embedding), collapse = ", ")))
      
      # Prepare plot data
      plot_data_runtime <- combined_runtime %>%
        mutate(
          n_sample = factor(n_sample, levels = n_samples),
          embedding = factor(embedding, levels = runtime_embeddings, 
                           labels = c("Scratch", "MedicalNet-ft"))
        )
      
      # Create color palette for runtime embeddings
      # Use colors 5 and 6 from the existing palette to match the other plots
      palet_discrete <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
      runtime_embedding_colors <- setNames(
        palet_discrete[c(6, 5)],  # Colors for Scratch and MedicalNet-ft
        c("Scratch", "MedicalNet-ft")
      )
      
      # Create boxplot
      p_runtime <- plot_data_runtime %>%
        ggplot(aes(x = n_sample, y = runtime, fill = embedding)) +
        geom_boxplot(position = position_dodge(width = 0.8), outlier.size = 0.8) +
        labs(
          x = "Sample Size",
          y = "Runtime (seconds)",
          fill = "Embedding"
        ) +
        theme_bw(base_size = 18) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 20),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          legend.position = "right",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16)
        ) +
        scale_fill_manual(values = runtime_embedding_colors) +
        scale_x_discrete(breaks = c("256", "460", "825", "1100", "5000")) +
        scale_y_log10(
          breaks = c(100, 500, 1000, 5000, 10000, 50000),
          labels = scales::comma
        )
      
      # Save runtime plots
      cat("\nSaving runtime boxplots:\n")
      
      # Save as PNG
      png_path_runtime <- file.path(output_dir, 
                                     sprintf("runtime_embedding_scratch_vs_medicalnet_ft_%d_%d.png", 
                                             min(seeds), max(seeds)))
      ggsave(
        png_path_runtime, 
        plot = p_runtime, 
        width = 12, 
        height = 7, 
        units = "in", 
        dpi = 300
      )
      cat("  Saved PNG:", png_path_runtime, "\n")
      
      # Save as PDF
      pdf_path_runtime <- file.path(output_dir, 
                                     sprintf("runtime_embedding_scratch_vs_medicalnet_ft_%d_%d.pdf", 
                                             min(seeds), max(seeds)))
      ggsave(
        pdf_path_runtime, 
        plot = p_runtime, 
        width = 12, 
        height = 7, 
        units = "in"
      )
      cat("  Saved PDF:", pdf_path_runtime, "\n")
      
      cat("\n=== Runtime Boxplots Complete ===\n")
      
    } else {
      warning("No runtime data loaded. Check if files exist in the runtime directory.")
    }
  }
} else {
  cat("\n=== Skipping Runtime Boxplots (only created for eps_sigmaY = 0.5) ===\n")
}

# ============================================================================
# Create QQ Plots for cVAE_long 200-seed Results (CI and No_CI)
# ============================================================================

cat("\n\n=== Creating QQ Plots for cVAE_long 200-seed Results ===\n")

# Load required libraries (should already be loaded, but ensure)
library(gridExtra)
library(grid)
library(cowplot)

# Define sample sizes for 200-seed analysis
n_samples_200 <- c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)

# Define subset of sample sizes for plotting (6 rows)
n_samples_plot_200 <- c(145, 256, 460, 1100, 5000, 10000)

# Map sample sizes to column indices (V1=145, V2=256, etc.)
sample_to_col_idx <- match(n_samples_plot_200, n_samples_200)

# Define column configurations: embedding-CIT-z combinations (6 columns)
column_configs <- list(
  list(cit = "RCOT_1", z_dim = "ukb_z1", label = "RCoT-z1"),
  list(cit = "RCOT_1", z_dim = "ukb_z6", label = "RCoT-z6"),
  list(cit = "RCOT_1", z_dim = "ukb_z15", label = "RCoT-z15"),
  list(cit = "comets_pcm", z_dim = "ukb_z1", label = "PCM-z1"),
  list(cit = "comets_pcm", z_dim = "ukb_z6", label = "PCM-z6"),
  list(cit = "comets_pcm", z_dim = "ukb_z15", label = "PCM-z15")
)

# Define conditions (CI and No_CI)
conditions_200 <- c("CI", "No_CI")

# Colors for the QQ plots
palet_qq <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
qq_color_condVAE <- palet_qq[3]  # condVAE color
#qq_color_cVAE_long <- palet_qq[7]  # cVAE_long color
qq_color_cVAE_long <- colorspace::lighten(palet_qq[7], amount = 0.2)

for (cond_200 in conditions_200) {
  cat(sprintf("\nCreating QQ plot grid for condition: %s\n", cond_200))
  
  # Determine the p-values path
  pval_path_200 <- file.path(results_base_path, cond_200, "p-values", "seeds_1_200")
  
  # Check if directory exists
  if (!dir.exists(pval_path_200)) {
    warning(sprintf("P-values directory not found: %s", pval_path_200))
    next
  }
  
  # Number of rows (sample sizes) and columns (CIT-z combinations)
  n_rows_qq <- length(n_samples_plot_200)
  n_cols_qq <- length(column_configs)
  
  # List to store all plots
  plot_list_200 <- list()
  
  # Loop over each sample size (rows)
  for (row_idx in seq_along(n_samples_plot_200)) {
    n_s <- n_samples_plot_200[row_idx]
    col_idx_data <- sample_to_col_idx[row_idx]  # Column index in data (1-10)
    
    # Loop over each column configuration
    for (col_idx in seq_along(column_configs)) {
      config <- column_configs[[col_idx]]
      cit_name <- config$cit
      z_dim <- config$z_dim
      col_label <- config$label
      
      # Construct filenames for both embeddings
      pval_filename_cVAE_long <- sprintf("1_0_1_0_fastsurfer_cVAE_long_%s_squared_%s.csv", 
                                          z_dim, cit_name)
      pval_filepath_cVAE_long <- file.path(pval_path_200, pval_filename_cVAE_long)
      
      pval_filename_condVAE <- sprintf("1_0_1_0_fastsurfer_condVAE_%s_squared_%s.csv", 
                                        z_dim, cit_name)
      pval_filepath_condVAE <- file.path(pval_path_200, pval_filename_condVAE)
      
      # Load cVAE_long p-values
      p_vals_cVAE_long <- NULL
      if (file.exists(pval_filepath_cVAE_long)) {
        pval_df <- fread(pval_filepath_cVAE_long, nThread = 1)
        pval_df <- pval_df[, -1]  # Drop first column
        col_name <- paste0("V", col_idx_data)
        p_vals_cVAE_long <- pval_df[[col_name]]
        p_vals_cVAE_long <- p_vals_cVAE_long[!is.na(p_vals_cVAE_long)]
      }
      
      # Load condVAE p-values
      p_vals_condVAE <- NULL
      if (file.exists(pval_filepath_condVAE)) {
        pval_df <- fread(pval_filepath_condVAE, nThread = 1)
        pval_df <- pval_df[, -1]  # Drop first column
        col_name <- paste0("V", col_idx_data)
        p_vals_condVAE <- pval_df[[col_name]]
        p_vals_condVAE <- p_vals_condVAE[!is.na(p_vals_condVAE)]
      }
      
      # Check if we have data for at least one embedding
      has_data <- (length(p_vals_cVAE_long) > 0 || length(p_vals_condVAE) > 0)
      
      if (has_data) {
        # Create combined data frame for plotting
        plot_data_list <- list()
        
        if (length(p_vals_cVAE_long) > 0) {
          ks_stat_cVAE_long <- round(ks.test(p_vals_cVAE_long, "punif")$statistic, 3)
          plot_data_list[[1]] <- data.frame(
            Observed = sort(p_vals_cVAE_long),
            Theoretical = qunif(ppoints(length(p_vals_cVAE_long))),
            Embedding = "cVAE_long"
          )
        }
        
        if (length(p_vals_condVAE) > 0) {
          ks_stat_condVAE <- round(ks.test(p_vals_condVAE, "punif")$statistic, 3)
          plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
            Observed = sort(p_vals_condVAE),
            Theoretical = qunif(ppoints(length(p_vals_condVAE))),
            Embedding = "condVAE"
          )
        }
        
        plot_df <- do.call(rbind, plot_data_list)
        
        # Construct KS statistic label
        if (length(p_vals_cVAE_long) > 0 && length(p_vals_condVAE) > 0) {
          ks_label <- sprintf("%.3f/%.3f", ks_stat_cVAE_long, ks_stat_condVAE)
        } else if (length(p_vals_cVAE_long) > 0) {
          ks_label <- sprintf("%.3f", ks_stat_cVAE_long)
        } else {
          ks_label <- sprintf("%.3f", ks_stat_condVAE)
        }
        
        # Create the QQ plot with both embeddings
        p <- ggplot() +
          geom_abline(slope = 1, intercept = 0, color = "black") +
          geom_point(data = plot_df, aes(x = Theoretical, y = Observed, color = Embedding), 
                     size = 1.5, alpha = 0.7) +
          scale_color_manual(values = c("cVAE_long" = qq_color_cVAE_long, 
                                        "condVAE" = qq_color_condVAE)) +
          xlab(ks_label) +
          ylab(n_s) +
          theme_minimal() +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
        
        # Base theme
        p <- p + theme(
          plot.title = element_blank(),
          text = element_text(size = 15),
          axis.title.y = element_blank(),
          legend.position = "none"  # Hide individual plot legends
        )
        
        # Visibility logic
        # First column: show y-axis with sample size label
        if (col_idx == 1) {
          p <- p + theme(
            axis.title.y = element_text(angle = 90, margin = margin(t = 0, r = 5, b = 0, l = 0))
          )
        } else {
          # Other columns: remove y-axis
          p <- p + theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
        }
        
        # Last row: show x-axis with KS statistic
        if (row_idx != n_rows_qq) {
          # Not last row: remove x-axis text and ticks
          p <- p + theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          )
        }
        
        plot_list_200[[length(plot_list_200) + 1]] <- p
      } else {
        # No data for either embedding - create empty plot
        cat(sprintf("  Warning: No data found for %s %s at n=%d\n", cit_name, z_dim, n_s))
        p <- ggplot() +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "No data", size = 5, hjust = 0.5)
        plot_list_200[[length(plot_list_200) + 1]] <- p
      }
    }
  }
  
  # Create column headers (CIT-z labels)
  col_headers <- lapply(column_configs, function(config) {
    textGrob(config$label, gp = gpar(fontsize = 15), just = "center")
  })
  
  # Create a legend grob
  # Create a dummy plot to extract legend
  dummy_data <- data.frame(
    x = c(1, 2),
    y = c(1, 2),
    Embedding = c("cVAE_long", "condVAE")
  )
  dummy_plot <- ggplot(dummy_data, aes(x = x, y = y, color = Embedding)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("cVAE_long" = qq_color_cVAE_long, 
                                   "condVAE" = qq_color_condVAE),
                       labels = c("cVAE_long" = "cVAE long", "condVAE" = "condVAE")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14))
  
  # Extract legend
  legend_grob <- cowplot::get_legend(dummy_plot)
  
  # Combine headers and plots
  grob_list <- c(col_headers, plot_list_200)
  
  # Create layout matrix:
  # First row: column headers
  # Subsequent rows: plots (sample sizes x column configs)
  layout_matrix <- rbind(
    matrix(1:n_cols_qq, nrow = 1, ncol = n_cols_qq),
    matrix((n_cols_qq + 1):(n_cols_qq + n_cols_qq * n_rows_qq), 
           nrow = n_rows_qq, ncol = n_cols_qq, byrow = TRUE)
  )
  
  # Arrange plots into a grid
  grid_qq_200_temp <- gridExtra::grid.arrange(
    grobs = grob_list,
    layout_matrix = layout_matrix,
    bottom = textGrob("Theoretical Quantiles", gp = gpar(fontsize = 20)),
    left = textGrob("Sample Quantiles", gp = gpar(fontsize = 20), rot = 90),
    heights = unit.c(unit(1, "lines"), unit(1.1, "null"), rep(unit(1, "null"), n_rows_qq - 1))
  )
  
  # Combine grid with legend at bottom
  grid_qq_200 <- gridExtra::grid.arrange(
    grid_qq_200_temp,
    legend_grob,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - unit(1.5, "lines"), unit(1.5, "lines"))
  )
  
  # Determine condition label for filename
  cond_label <- ifelse(cond_200 == "CI", "CI", "No_CI")
  
  # Save QQ plot
  cat(sprintf("Saving QQ plot for %s:\n", cond_label))
  
  # Save as PNG
  png_path_qq_200 <- file.path(output_dir, 
                                sprintf("qq_plot_cVAE_long_%s_seeds_1_200.png", cond_label))
  ggsave(
    png_path_qq_200, 
    plot = grid_qq_200, 
    width = n_cols_qq * 2.5, 
    height = n_rows_qq * 2.5, 
    units = "in", 
    dpi = 300
  )
  cat("  Saved PNG:", png_path_qq_200, "\n")
  
  # Save as PDF
  pdf_path_qq_200 <- file.path(output_dir, 
                                sprintf("qq_plot_cVAE_long_%s_seeds_1_200.pdf", cond_label))
  ggsave(
    pdf_path_qq_200, 
    plot = grid_qq_200, 
    width = n_cols_qq * 2.5, 
    height = n_rows_qq * 2.5, 
    units = "in"
  )
  cat("  Saved PDF:", pdf_path_qq_200, "\n")
}

cat("\n=== QQ Plots for cVAE_long 200-seed Results Complete ===\n")
cat(sprintf("Output directory: %s\n", output_dir))

### Runtime condVAE vs cVAE_long
results_base_path <- "/sc/home/marco.simnacher/dncitPaper/Results"
# Determine the runtime path
runtime_path_200 <- file.path(results_base_path, "CI", "runtime_cit", "seeds_1_200")
# Construct filenames for both embeddings
runtime_filename_cVAE_long_rcot <- sprintf("1_0_1_0_fastsurfer_cVAE_long_%s_squared_%s.csv", 
                                          "ukb_z1", "RCOT_1")
runtime_filepath_cVAE_long_rcot <- file.path(runtime_path_200, runtime_filename_cVAE_long_rcot)
      
runtime_filename_condVAE_rcot <- sprintf("1_0_1_0_fastsurfer_condVAE_%s_squared_%s.csv", 
                                        "ukb_z1", "RCOT_1")
runtime_filepath_condVAE_rcot <- file.path(runtime_path_200, runtime_filename_condVAE_rcot)

# Load cVAE_long runtime
runtime_df_cVAE_long_rcot <- fread(runtime_filepath_cVAE_long_rcot, nThread = 1)
runtime_df_cVAE_long_rcot <- runtime_df_cVAE_long_rcot[, -1]
runtime_df_condVAE_rcot <- fread(runtime_filepath_condVAE_rcot, nThread = 1)
runtime_df_condVAE_rcot <- runtime_df_condVAE_rcot[, -1]
rcot_cVAE_long_mean <- colMeans(runtime_df_cVAE_long_rcot)
rcot_condVAE_mean <- colMeans(runtime_df_condVAE_rcot)
rcot_runtime_mean <- rcot_cVAE_long_mean - rcot_condVAE_mean
rcot_runtime_mean


runtime_filename_cVAE_long_pcm <- sprintf("1_0_1_0_fastsurfer_cVAE_long_%s_squared_%s.csv", 
                                          "ukb_z1", "comets_pcm")
runtime_filepath_cVAE_long_pcm <- file.path(runtime_path_200, runtime_filename_cVAE_long_pcm)
      
runtime_filename_condVAE_pcm <- sprintf("1_0_1_0_fastsurfer_condVAE_%s_squared_%s.csv", 
                                        "ukb_z1", "comets_pcm")
runtime_filepath_condVAE_pcm <- file.path(runtime_path_200, runtime_filename_condVAE_pcm)

# Load cVAE_long runtime
runtime_df_cVAE_long_pcm <- fread(runtime_filepath_cVAE_long_pcm, nThread = 1)
runtime_df_cVAE_long_pcm <- runtime_df_cVAE_long_pcm[, -1]
runtime_df_condVAE_pcm <- fread(runtime_filepath_condVAE_pcm, nThread = 1)
runtime_df_condVAE_pcm <- runtime_df_condVAE_pcm[, -1]
pcm_cVAE_long_mean <- colMeans(runtime_df_cVAE_long_pcm)
pcm_condVAE_mean <- colMeans(runtime_df_condVAE_pcm)
pcm_runtime_mean <- pcm_cVAE_long_mean - pcm_condVAE_mean
pcm_runtime_mean
sample_sizes <- c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)
cor(rcot_runtime_mean, sample_sizes)
cor(pcm_runtime_mean, sample_sizes)
#scatter plot runtime_mean vs sample sizes
pdf_path_runtime <- file.path(output_dir, "runtime_comparison_rcot_pcm.pdf")
pdf(pdf_path_runtime, width = 10, height = 7)
plot(sample_sizes, pcm_runtime_mean, type = "l", col = "red", xlab = "Sample Size", ylab = "Runtime (seconds)")
lines(sample_sizes, rcot_runtime_mean, type = "l", col = "blue")
legend("topright", legend = c("RCoT", "PCM"), col = c("red", "blue"), lty = 1)
dev.off()
