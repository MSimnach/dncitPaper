library(dncitPaper)
dncits <- c("RCOT", "CMIknn", "kpc_graph_2_10", "FCIT", "WALD" )
settings <- data.frame(
    confounder = rep(c("age, sex", paste("age, sex, 10 PCs")), each=12),
    fct_relation = rep(rep(c("linear", "non_lin", "strong_non_lin"), each=4), 2),
    embedding_map = rep(c("same", "noisy", "noisy", "different"), 6),
    noise_x = rep(c(0, 50,500, 0), 6)
)

#Type 1 rejection rates
folder_path_T1E <- "C:/Users/Marco/seadrive_root/Marco Si/Meine Bibliotheken/Meine Bibliothek/CITs/P1 Deep CITs/Paper_DNCITs/Results/CI/rejection_rates"
sample_sizes <- c(30, 100, 300, 1000, 3000, 10000)
p_T1E <- create_test_plot(folder_path_T1E, dncits, settings, sample_sizes = sample_sizes)

# Power rejection rates
folder_path_power <-  "C:/Users/Marco/seadrive_root/Marco Si/Meine Bibliotheken/Meine Bibliothek/CITs/P1 Deep CITs/Paper_DNCITs/Results/No_CI/rejection_rates"
betas <- c(100.00, 75.00, 50.00, 25.00, 10.00, 7.50, 5.00, 2.50, 1.00, 0.75, 0.50, 0.25)
p_power <- create_test_plot(folder_path_power, dncits, settings, betas=betas)

#Avg runtimes for each setting
folder_path_runtime <- "C:/Users/Marco/seadrive_root/Marco Si/Meine Bibliotheken/Meine Bibliothek/CITs/P1 Deep CITs/Paper_DNCITs/Results/CI/Runtime"
sample_sizes <- c(30, 100, 300, 1000, 3000, 10000)
runtimes <- create_test_plot(folder_path_runtime, dncits, settings, sample_sizes = sample_sizes, runtimes_plot=TRUE)



#Avg runtimes averaged over all settings
##color palettes for plotting
palette_discrete <- paletteer::paletteer_c("grDevices::RdYlBu", 20)
folder_path_runtime <- "C:/Users/Marco/seadrive_root/Marco Si/Meine Bibliotheken/Meine Bibliothek/CITs/P1 Deep CITs/Paper_DNCITs/Results/CI/Runtime"

all_files <- list.files(folder_path_runtime, full.names = TRUE)

settings_files <- unique(substr(basename(all_files), 1, 39))
# exclude additional noise for Y
settings_files <- settings_files[!grepl('_50_0_fas', settings_files)]
# exclude noise 10 on embeddings, since tests not sensitive for noise on features
settings_files <- settings_files[!grepl('10_1_0_fas', settings_files)]
#Reorder matching with tree diagram in y axis
settings_files <- settings_files[c(23,21,17,19,15,13,9,11,7,5,1,3,24,22,18,20,16,14,10,12,8,6,2,4)]

# check if all DNCITs have results for all settings
for(dncit in dncits){
  files_with_CMIknn <- list.files(path = folder_path_runtime, pattern = dncit, full.names = TRUE)
  print(paste0(dncit, " has results in ", length(files_with_CMIknn), " settings."))
}

# empty result tab for each DNCIT
list_result_tabs <- list()
betas <- NULL
for(dncit in dncits){
  dncit_char <- paste0(dncit)
  list_result_tabs[[dncit_char]] <- create_empty_result_tab(settings_files,sample_sizes=sample_sizes, betas=betas)
}

list_result_tabs_full <- create_runtime_tabs(dncits, settings_files, all_files, list_result_tabs, sample_sizes)

avg_runtime <- list()
for(dncit in dncits){
  tmp <- list_result_tabs_full[[dncit]][1:12,2:7]
  tmp[] <- lapply(list_result_tabs_full[[dncit]][1:12,2:7], function(x) as.numeric(as.character(x)))
  avg_runtime[[paste0(dncit, '_AS')]] <- colMeans(tmp, na.rm=TRUE)

  tmp <- list_result_tabs_full[[dncit]][13:24,2:7]
  tmp[] <- lapply(list_result_tabs_full[[dncit]][13:24,2:7], function(x) as.numeric(as.character(x)))
  avg_runtime[[paste0(dncit, '_genes')]] <- colMeans(tmp, na.rm=TRUE)
}
runtime_tbl <- tibble::as_tibble(avg_runtime) %>%
  mutate(sample_sizes = sample_sizes)

runtime_long <- runtime_tbl %>%
  tidyr::pivot_longer(cols = -sample_sizes, names_to = "measurement", values_to = "value")

new_names <- c(RCOT_AS = "rcot 2", CMIknn_AS = "CMIknn 2", kpc_graph_2_10_AS = "cpt-kpc 2", FCIT_AS = "fcit 2", WALD_AS = "Wald 2",
               RCOT_genes = "rcot 12", CMIknn_genes = "CMIknn 12", kpc_graph_2_10_genes = "cpt-kpc 12", FCIT_genes = "fcit 12", WALD_genes = "Wald 12")
runtime_long <- runtime_long %>%
  mutate(measurement = recode(measurement, !!!new_names))

# View reshaped data
ggplot2::ggplot(runtime_long, ggplot2::aes(x = log(sample_sizes), y = value, color = measurement,group = measurement)) +
  ggplot2::geom_line() +  # Add line to connect points
  ggplot2::geom_point() +  # Add points on each data point
  ggplot2::scale_color_manual(values = palette_discrete[c(1,2,4,5,7,8,15,16,18,20)]) +  # Apply the specified color palette
  ggplot2::scale_x_continuous(
    name = "Sample Sizes",
    breaks = log(runtime_tbl$sample_sizes),
    labels = runtime_tbl$sample_sizes
  ) +
  ggplot2::labs(x = "Sample Sizes", y = "Runtime") +
  ggplot2::theme_minimal() +
  ggplot2::theme()
ggplot2::ggsave("C:/Users/Marco/seadrive_root/Marco Si/Meine Bibliotheken/Meine Bibliothek/CITs/P1 Deep CITs/Paper_DNCITs/plots_simulation/runtime.png", width = 7, height = 7, dpi = 300)


