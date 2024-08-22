##### Plot rejection rates of CITs in nested loop plots (cf. Rücker et al. 2014)
## For presentations
##color palettes for plotting
palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")
path_to_save_nested_loop_plots <- "C:\\Users\\Marco\\seadrive_root\\Marco Si\\Meine Bibliotheken\\Meine Bibliothek\\CITs\\P1 Deep CITs\\Paper_DNCITs\\Plots_simulation\\"

#### split into CI and No CI
folder_path <- "M:\\CITs\\Application\\UKB_data\\Results\\No_CI\\rejection_rates"
#folder_path <- "M:\\CITs\\Application\\UKB_data\\Results\\CI\\rejection_rates"
all_files <- list.files(folder_path, full.names = TRUE)
all_files <- all_files[setdiff(1:length(all_files), grep('2_0_1_0|2_3_1_0|3_0_1_0|3_3_1_0|4_0_1_0|4_3_1_0|5_0_1_0|5_3_1_0', all_files))]

### 1) Data preparation for nested loop over
# loop 1: confounder dimension (1,2,4,6,10,15)
# loop 2: sample size
# fixed confounder relationship (squared terms of all continuous confounders)
# Deep-WALD and Deep-RCOT for all four embedding maps (freesurfer, fastsurfer, condVAE, noisy fastsurfer)
squared_conf_files <- all_files[grep("squared", all_files)]
wald_files <- all_files[grep("WALD", all_files)]
rcot_files <- all_files[grep("RCOT", all_files)]
files_conf_dim <- union(intersect(squared_conf_files, wald_files), intersect(squared_conf_files, rcot_files))

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
      confounder = c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15'))
design <- expand.grid(params)
design <- design %>%
  mutate("fastsurfer_fastsurfer-RCOT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-WALD" = rep(NA, nrow(design)),
         "freesurfer-RCOT" = rep(NA, nrow(design)),
         "freesurfer-WALD" = rep(NA, nrow(design)),
         "condVAE-RCOT" = rep(NA, nrow(design)),
         "condVAE-WALD" = rep(NA, nrow(design)))

embedding_maps <- c('fastsurfer_fastsurfer', 'freesurfer', 'condVAE')
dncits <- c('RCOT', 'WALD')
confounder <- c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15')
for(dncit in dncits){
  for (embedding in embedding_maps){
    for (conf in confounder){
      files_dncit <- grep(dncit, files_conf_dim, value = TRUE)
      files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
      files <- grep(conf, files_dncit_conf, value=TRUE)
      df <- read.csv(files, header = TRUE, sep = ",")
      col <- grepl(embedding, colnames(design)) & grepl(dncit, colnames(design))
      design[design$confounder==conf, col] <- df[,2]
    }
  }
}

design$confounder <- rep(c(1,2,4,6,10,15), each=8)
colnames(design) <- c("sample_sizes", "confounder dimension", "Fastsurfer-RCOT", "Fastsurfer-WALD",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "condVAE-RCOT", "condVAE-WALD")

design_conf_dim_ci <- design
#design_conf_dim_no_ci <- design

### 2) Data preparation for nested loop over
# loop 1: confounder functional relationship g_z ("linear", "squared", "realistic")
# loop 2: sample size
# fixed confounder dimension 6
# Deep-WALD and Deep-RCOT for all four embedding maps (freesurfer, fastsurfer, condVAE, noisy fastsurfer)
ukb_z6_conf_files <- all_files[grep("ukb_z6", all_files)]
wald_files <- all_files[grep("WALD", all_files)]
rcot_files <- all_files[grep("RCOT", all_files)]
files_conf_relation <- union(intersect(ukb_z6_conf_files, wald_files), intersect(ukb_z6_conf_files, rcot_files))

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
               confounder = c('linear', 'squared', 'realistic'))#, 'breakpoint3'))
design <- expand.grid(params)
design <- design %>%
  mutate("fastsurfer_fastsurfer-RCOT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-WALD" = rep(NA, nrow(design)),
         "freesurfer-RCOT" = rep(NA, nrow(design)),
         "freesurfer-WALD" = rep(NA, nrow(design)),
         "condVAE-RCOT" = rep(NA, nrow(design)),
         "condVAE-WALD" = rep(NA, nrow(design)))

embedding_maps <- c('fastsurfer_fastsurfer', 'freesurfer', 'condVAE')
dncits <- c('RCOT', 'WALD')
confounder <- c('linear', 'squared', 'realistic')#, 'breakpoint3')
for(dncit in dncits){
  for (embedding in embedding_maps){
    for (conf in confounder){
      files_dncit <- grep(dncit, files_conf_relation, value = TRUE)
      files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
      files <- grep(conf, files_dncit_conf, value=TRUE)
      df <- read.csv(files, header = TRUE, sep = ",")
      col <- grepl(embedding, colnames(design)) & grepl(dncit, colnames(design))
      design[design$confounder==conf, col] <- df[,2]
    }
  }
}

design$confounder <- rep(c('linear', 'squared', 'realistic'), each=8)
colnames(design) <- c("sample_sizes", "confounder", "Fastsurfer-RCOT", "Fastsurfer-WALD",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "condVAE-RCOT", "condVAE-WALD")

#design_conf_relation_ci <- design
design_conf_relation_no_ci <- design



##### nested loop plots
## 1) splitted into T1E and power
## CI
p_ci <- looplot::nested_loop_plot(resdf = design,
                              x = "sample_sizes", steps = "confounder dimension",
                              steps_y_base = -0.05, steps_y_height = 0.05,
                              x_name = "Sample size", y_name = "Rejection rate",
                              spu_x_shift = 1,
                              colors = palet_discrete[rep(c(2,9),3)],
                              line_linetypes = rep(c(1,2,3), each=2),
                              sizes = rep(2,6),
                              line_size = 1.5,
                              point_shapes = rep(c(15,19),3),
                              steps_values_annotate = TRUE, steps_annotation_size = 4,
                              hline_intercept = c(0, 0.05),
                              hline_linetype = c(3,1),
                              hline_size = c(1,1.5),
                              y_expand_add = c(0.05,0.05),
                              legend_name = "DNCIT",
                              post_processing = list(
                                add_custom_theme = list(
                                  axis.text.x = ggplot2::element_text(angle = -90,
                                                                      vjust = 0.5,
                                                                      size = 8)
                                )
                              ))
print(p_ci)

## No CI
p_no_ci = looplot::nested_loop_plot(resdf = design,
                              x = "sample_sizes", steps = "confounder dimension",
                              steps_y_base = -0.05, steps_y_height = 0.05,
                              x_name = "Sample size", y_name = "Rejection rate",
                              spu_x_shift = 1,
                              colors = palet_discrete[rep(c(2,9),3)],
                              line_linetypes = rep(c(1,2,3), each=2),
                              sizes = rep(2,6),
                              line_size = 1.5,
                              point_shapes = rep(c(15,19),3),
                              steps_values_annotate = TRUE, steps_annotation_size = 4,
                              hline_intercept = 0,
                              y_expand_add = c(0.05,0.05),
                              legend_name = "DNCIT",
                              post_processing = list(
                                add_custom_theme = list(
                                  axis.text.x = ggplot2::element_text(angle = -90,
                                                                      vjust = 0.5,
                                                                      size = 8)
                                )
                              ))
print(p_no_ci)

## 2) power and T1E together
# conf relation
design <- rbind(design_conf_relation_ci, design_conf_relation_no_ci)
design$'Setting' <- rep(c('No', 'Yes'), each=24)
design <- design %>%
  mutate(confounder = recode(confounder,
                             'linear' = '1',
                             'squared' = '2',
                             'realistic' = '3'))#,
                             #'complex' = 4))
p_conf_relation <- looplot::nested_loop_plot(resdf = design,
                                    x = "sample_sizes",
                                    grid_rows = 'Setting',
                                    steps = "confounder",
                                    steps_y_base = -0.1, steps_y_height = 0.05,
                                    x_name = "Sample size", y_name = "Rejection rate",
                                    spu_x_shift = 1,
                                    colors = palet_discrete[rep(c(2,9),3)],
                                    line_linetypes = rep(c(1,3,6), each=2),
                                    line_size = 1.5,
                                    point_shapes = rep(c(15,19),3),
                                    steps_values_annotate = TRUE, steps_annotation_size = 6,
                                    hline_intercept = 0,
                                    y_expand_add = c(0.1,0.15),
                                    legend_name = "DNCIT",
                                    base_size = 24,
                                    replace_labels = list(
                                      Setting = c('No'='T1E',
                                                     'Yes'='Power'),
                                      confounder = c('1'='linear',
                                                     '2'='squared',
                                                     '3'='complex')
                                    ),
                                    post_processing = list(
                                      add_custom_theme = list(
                                        axis.text.x = ggplot2::element_text(angle = -90,
                                                                            vjust = 0.5,
                                                                            size = 15)
                                      )
                                    ))
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, 'nested_loop_conf_rel_power_T1E.png'), p_conf_relation, width = 16, height = 10, dpi = 300)
print(p_conf_relation)

# conf dimension
design <- rbind(design_conf_dim_ci, design_conf_dim_no_ci)
design$'Setting' <- rep(c('No', 'Yes'), each=48)
p_conf_dim = looplot::nested_loop_plot(resdf = design,
                                            x = "sample_sizes",
                                            grid_rows = 'Setting',
                                            steps  = "confounder dimension",
                                            steps_y_base = -0.1, steps_y_height = 0.05,
                                            x_name = "Sample size", y_name = "Rejection rate",
                                            spu_x_shift = 1,
                                            colors = palet_discrete[rep(c(2,9),3)],
                                            line_linetypes = rep(c(1,3,6), each=2),
                                            #sizes = rep(3,6),
                                            line_size = 1.5,
                                            point_shapes = rep(c(15,19),3),
                                            steps_values_annotate = TRUE, steps_annotation_size = 6,
                                            hline_intercept = 0,
                                            y_expand_add = c(0.1,0.15),
                                            legend_name = "DNCIT",
                                       base_size = 24,
                                            replace_labels = list(
                                              Setting = c('No'='T1E',
                                                            'Yes'='Power')
                                            ),
                                            post_processing = list(
                                              add_custom_theme = list(
                                                axis.text.x = ggplot2::element_text(angle = -90,
                                                                                    vjust = 0.5,
                                                                                    size = 15)
                                              )
                                            ))
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, 'nested_loop_conf_dim.png'), p_conf_dim, width = 16, height = 10, dpi = 300)
print(p_conf_dim)


## For the paper
##color palettes for plotting
palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")
path_to_save_nested_loop_plots <- "C:\\Users\\Marco\\seadrive_root\\Marco Si\\Meine Bibliotheken\\Meine Bibliothek\\CITs\\P1 Deep CITs\\Paper_DNCITs\\Plots_simulation\\"

#### split into CI and No CI
folder_path <- "M:\\CITs\\Application\\UKB_data\\Results\\No_CI\\rejection_rates"
#folder_path <- "M:\\CITs\\Application\\UKB_data\\Results\\CI\\rejection_rates"
all_files <- list.files(folder_path, full.names = TRUE)
all_files <- all_files[setdiff(1:length(all_files), grep('2_0_1_0|2_3_1_0|3_0_1_0|3_3_1_0|4_0_1_0|4_3_1_0|5_0_1_0|5_3_1_0', all_files))]

### 1) Data preparation for nested loop over
# loop 1: confounder dimension (1,2,4,6,10,15)
# loop 2: sample size
# fixed confounder relationship (squared terms of all continuous confounders)
# Deep-WALD and Deep-RCOT for all four embedding maps (freesurfer, fastsurfer, condVAE, noisy fastsurfer)
squared_conf_files <- all_files[grep("squared", all_files)]
wald_files <- all_files[grep("WALD", all_files)]
rcot_files <- all_files[grep("RCOT", all_files)]
kpc_cpt_files <- all_files[grep("kpc_graph", all_files)]
files_conf_dim <- c(intersect(squared_conf_files, wald_files), intersect(squared_conf_files, rcot_files),
                        intersect(squared_conf_files, kpc_cpt_files))

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
               confounder = c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15'))
design <- expand.grid(params)
design <- design %>%
  mutate("fastsurfer_fastsurfer-RCOT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-WALD" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-kpc_graph" = rep(NA, nrow(design)),
         "freesurfer-RCOT" = rep(NA, nrow(design)),
         "freesurfer-WALD" = rep(NA, nrow(design)),
         "freesurfer-kpc_graph" = rep(NA, nrow(design)),
         "condVAE-RCOT" = rep(NA, nrow(design)),
         "condVAE-WALD" = rep(NA, nrow(design)),
         "condVAE-kpc_graph" = rep(NA, nrow(design)))

embedding_maps <- c('fastsurfer_fastsurfer', 'freesurfer', 'condVAE')
dncits <- c('RCOT', 'WALD', 'kpc_graph')
confounder <- c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15')
for(dncit in dncits){
  for (embedding in embedding_maps){
    for (conf in confounder){
      files_dncit <- grep(dncit, files_conf_dim, value = TRUE)
      files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
      files <- grep(conf, files_dncit_conf, value=TRUE)
      df <- read.csv(files, header = TRUE, sep = ",")
      col <- grepl(embedding, colnames(design)) & grepl(dncit, colnames(design))
      design[design$confounder==conf & design$sample_sizes %in% df[,1], col] <- df[,2]
    }
  }
}

design$confounder <- rep(c(1,2,4,6,10,15), each=10)
colnames(design) <- c("sample_sizes", "confounder dimension", "Fastsurfer-RCOT", "Fastsurfer-WALD", "Fastsurfer-CPT_KPC",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "Freesurfer-CPT_KPC", "condVAE-RCOT", "condVAE-WALD", "condVAE-CPT_KPC")

#design_conf_dim_ci <- design
design_conf_dim_no_ci <- design

### 2) Data preparation for nested loop over
# loop 1: confounder functional relationship g_z ("linear", "squared", "realistic")
# loop 2: sample size
# fixed confounder dimension 6
# Deep-WALD and Deep-RCOT for all four embedding maps (freesurfer, fastsurfer, condVAE, noisy fastsurfer)
ukb_z6_conf_files <- all_files[grep("ukb_z6", all_files)]
wald_files <- all_files[grep("WALD", all_files)]
rcot_files <- all_files[grep("RCOT", all_files)]
kpc_cpt_files <- all_files[grep("kpc_graph", all_files)]
files_conf_relation <- c(intersect(ukb_z6_conf_files, wald_files), intersect(ukb_z6_conf_files, rcot_files),
                         intersect(ukb_z6_conf_files, kpc_cpt_files))

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
               confounder = c('linear', 'squared', 'realistic'))#, 'breakpoint3'))
design <- expand.grid(params)
design <- design %>%
  mutate("fastsurfer_fastsurfer-RCOT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-WALD" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-kpc_graph" = rep(NA, nrow(design)),
         "freesurfer-RCOT" = rep(NA, nrow(design)),
         "freesurfer-WALD" = rep(NA, nrow(design)),
         "freesurfer-kpc_graph" = rep(NA, nrow(design)),
         "condVAE-RCOT" = rep(NA, nrow(design)),
         "condVAE-WALD" = rep(NA, nrow(design)),
         "condVAE-kpc_graph" = rep(NA, nrow(design)))

embedding_maps <- c('fastsurfer_fastsurfer', 'freesurfer', 'condVAE')
dncits <- c('RCOT', 'WALD', 'kpc_graph')
confounder <- c('linear', 'squared', 'realistic')#, 'breakpoint3')
for(dncit in dncits){
  for (embedding in embedding_maps){
    for (conf in confounder){
      files_dncit <- grep(dncit, files_conf_relation, value = TRUE)
      files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
      files <- grep(conf, files_dncit_conf, value=TRUE)
      df <- read.csv(files, header = TRUE, sep = ",")
      col <- grepl(embedding, colnames(design)) & grepl(dncit, colnames(design))
      design[design$confounder==conf& design$sample_sizes %in% df[,1], col] <- df[,2]
    }
  }
}

design$confounder <- rep(c('linear', 'squared', 'realistic'), each=10)
colnames(design) <- c("sample_sizes", "confounder", "Fastsurfer-RCOT", "Fastsurfer-WALD", "Fastsurfer-CPT_KPC",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "Freesurfer-CPT_KPC", "condVAE-RCOT", "condVAE-WALD", "condVAE-CPT_KPC")

#design_conf_relation_ci <- design
design_conf_relation_no_ci <- design



##### nested loop plots
## 1) splitted into T1E and power
## CI
p_ci <- looplot::nested_loop_plot(resdf = design,
                                  x = "sample_sizes", steps = "confounder dimension",
                                  steps_y_base = -0.05, steps_y_height = 0.05,
                                  x_name = "Sample size", y_name = "Rejection rate",
                                  spu_x_shift = 1,
                                  colors = palet_discrete[rep(c(2,9),3)],
                                  line_linetypes = rep(c(1,2,3), each=2),
                                  sizes = rep(2,6),
                                  line_size = 1.5,
                                  point_shapes = rep(c(15,19),3),
                                  steps_values_annotate = TRUE, steps_annotation_size = 4,
                                  hline_intercept = c(0, 0.05),
                                  hline_linetype = c(3,1),
                                  hline_size = c(1,1.5),
                                  y_expand_add = c(0.05,0.05),
                                  legend_name = "DNCIT",
                                  post_processing = list(
                                    add_custom_theme = list(
                                      axis.text.x = ggplot2::element_text(angle = -90,
                                                                          vjust = 0.5,
                                                                          size = 8)
                                    )
                                  ))
print(p_ci)

## No CI
p_no_ci = looplot::nested_loop_plot(resdf = design,
                                    x = "sample_sizes", steps = "confounder dimension",
                                    steps_y_base = -0.05, steps_y_height = 0.05,
                                    x_name = "Sample size", y_name = "Rejection rate",
                                    spu_x_shift = 1,
                                    colors = palet_discrete[rep(c(2,9),3)],
                                    line_linetypes = rep(c(1,2,3), each=2),
                                    sizes = rep(2,6),
                                    line_size = 1.5,
                                    point_shapes = rep(c(15,19),3),
                                    steps_values_annotate = TRUE, steps_annotation_size = 4,
                                    hline_intercept = 0,
                                    y_expand_add = c(0.05,0.05),
                                    legend_name = "DNCIT",
                                    post_processing = list(
                                      add_custom_theme = list(
                                        axis.text.x = ggplot2::element_text(angle = -90,
                                                                            vjust = 0.5,
                                                                            size = 8)
                                      )
                                    ))
print(p_no_ci)

## 2) power and T1E together
# conf relation
design <- rbind(design_conf_relation_ci, design_conf_relation_no_ci)
design$'Setting' <- rep(c('No', 'Yes'), each=30)
design <- design %>%
  mutate(confounder = recode(confounder,
                             'linear' = '1',
                             'squared' = '2',
                             'realistic' = '3'))#,
#'complex' = 4))
p_conf_relation <- looplot::nested_loop_plot(resdf = design,
                                             x = "sample_sizes",
                                             grid_rows = 'Setting',
                                             steps = "confounder",
                                             steps_y_base = -0.1, steps_y_height = 0.05,
                                             x_name = "Sample size", y_name = "Rejection rate",
                                             spu_x_shift = 1,
                                             colors = palet_discrete[rep(c(2, 4, 9),3)],
                                             line_linetypes = rep(c(1,3,6), each=3),
                                             line_size = 1.5,
                                             point_shapes = rep(c(15,17,19),3),
                                             steps_values_annotate = TRUE, steps_annotation_size = 6,
                                             hline_intercept = 0,
                                             y_expand_add = c(0.1,0.15),
                                             legend_name = "DNCIT",
                                             base_size = 24,
                                             replace_labels = list(
                                               Setting = c('No'='T1E',
                                                           'Yes'='Power'),
                                               confounder = c('1'='linear',
                                                              '2'='squared',
                                                              '3'='complex')
                                             ),
                                             post_processing = list(
                                               add_custom_theme = list(
                                                 axis.text.x = ggplot2::element_text(angle = -90,
                                                                                     vjust = 0.5,
                                                                                     size = 15)
                                               )
                                             ))
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, 'all_dncits_nested_loop_conf_rel_power_T1E.png'), p_conf_relation, width = 16, height = 10, dpi = 300)
print(p_conf_relation)

# conf dimension
design <- rbind(design_conf_dim_ci, design_conf_dim_no_ci)
design$'Setting' <- rep(c('No', 'Yes'), each=60)
p_conf_dim = looplot::nested_loop_plot(resdf = design,
                                       x = "sample_sizes",
                                       grid_rows = 'Setting',
                                       steps  = "confounder dimension",
                                       steps_y_base = -0.1, steps_y_height = 0.05,
                                       x_name = "Sample size", y_name = "Rejection rate",
                                       spu_x_shift = 1,
                                       colors = palet_discrete[rep(c(2,4,9),3)],
                                       line_linetypes = rep(c(1,3,6), each=3),
                                       #sizes = rep(3,6),
                                       line_size = 1.5,
                                       point_shapes = rep(c(15,17,19),3),
                                       steps_values_annotate = TRUE, steps_annotation_size = 6,
                                       hline_intercept = 0,
                                       y_expand_add = c(0.1,0.15),
                                       legend_name = "DNCIT",
                                       base_size = 24,
                                       replace_labels = list(
                                         Setting = c('No'='T1E',
                                                     'Yes'='Power')
                                       ),
                                       post_processing = list(
                                         add_custom_theme = list(
                                           axis.text.x = ggplot2::element_text(angle = -90,
                                                                               vjust = 0.5,
                                                                               size = 15)
                                         )
                                       ))
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, 'all_dncits_nested_loop_conf_dim.png'), p_conf_dim, width = 16, height = 10, dpi = 300)
print(p_conf_dim)

