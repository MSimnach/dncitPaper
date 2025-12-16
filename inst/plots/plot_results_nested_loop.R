##### Plot rejection rates of CITs in nested loop plots (cf. Rücker et al. 2014)
library(dplyr)
library(looplot)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(grid)
#### For presentations
##color palettes for plotting
palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")
path_to_save_nested_loop_plots <- "/sc/home/marco.simnacher/dncitPaper/inst/diagnostic_embedding/figures"

#### split into CI and No CI
#path to results (no ci and ci rejection rates)
folder_path <- "Results\\No_CI\\rejection_rates"
#folder_path <- "Results\\CI\\rejection_rates"
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
                                    base_size = 34,
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
                                                                            size = 22)
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
                                       base_size = 34,
                                            replace_labels = list(
                                              Setting = c('No'='T1E',
                                                            'Yes'='Power')
                                            ),
                                            post_processing = list(
                                              add_custom_theme = list(
                                                axis.text.x = ggplot2::element_text(angle = -90,
                                                                                    vjust = 0.5,
                                                                                    size = 14)
                                              )
                                            ))
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, 'nested_loop_conf_dim.png'), p_conf_dim, width = 16, height = 10, dpi = 600)
print(p_conf_dim)




#### For the paper
##color palettes for plotting
#palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")
palet_discrete <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
#palet_discrete <- palet_discrete[seq(1,15, length.out=12)]
path_to_save_nested_loop_plots <- "/sc/home/marco.simnacher/dncitPaper/inst/plots"

#### split into CI and No CI
folder_path <- "/sc/home/marco.simnacher/dncitPaper/Results/No_CI/rejection_rates/seeds_1_200"
#folder_path <- "Results\\No_CI\\rejection_rates"
all_files <- list.files(folder_path, full.names = TRUE)
all_files <- all_files[setdiff(1:length(all_files), grep('1_0_0.5_0', all_files))]

### 1) Data preparation for nested loop over
# loop 1: confounder dimension (1,2,4,6,10,15) [APPENDIX] vs confounder dimension (1,2,10) [MAIN TEXT]
# loop 2: sample size
# fixed confounder relationship (squared terms of all continuous confounders)
cit_patterns <- "WALD|RCOT_1.csv|kpc_graph|FCIT|CMIknn|comets_pcm"
cit_files <- grep(cit_patterns, all_files, value=TRUE)
files_conf_dim <- grep("squared", cit_files, value=TRUE)

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
               confounder = c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15'))
design <- expand.grid(params)
design <- design %>%
  mutate("fastsurfer_fastsurfer-RCOT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-WALD" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-kpc_graph" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-FCIT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-CMIknn" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-comets_pcm" = rep(NA, nrow(design)),
         "freesurfer-RCOT" = rep(NA, nrow(design)),
         "freesurfer-WALD" = rep(NA, nrow(design)),
         "freesurfer-kpc_graph" = rep(NA, nrow(design)),
         "freesurfer-FCIT" = rep(NA, nrow(design)),
         "freesurfer-CMIknn" = rep(NA, nrow(design)),
         "freesurfer-comets_pcm" = rep(NA, nrow(design)),
         "condVAE-RCOT" = rep(NA, nrow(design)),
         "condVAE-WALD" = rep(NA, nrow(design)),
         "condVAE-kpc_graph" = rep(NA, nrow(design)),
         "condVAE-FCIT" = rep(NA, nrow(design)),
         "condVAE-CMIknn" = rep(NA, nrow(design)),
         "condVAE-comets_pcm" = rep(NA, nrow(design)),
         "medicalnet-RCOT" = rep(NA, nrow(design)),
         "medicalnet-WALD" = rep(NA, nrow(design)),
         "medicalnet-kpc_graph" = rep(NA, nrow(design)),
         "medicalnet-FCIT" = rep(NA, nrow(design)),
         "medicalnet-CMIknn" = rep(NA, nrow(design)),
         "medicalnet-comets_pcm" = rep(NA, nrow(design)))

embedding_maps <- c('fastsurfer_fastsurfer', 'freesurfer', 'condVAE', 'medicalnet')
dncits <- c('RCOT', 'WALD', 'kpc_graph', 'FCIT', 'CMIknn', 'comets_pcm')
# [APPENDIX]
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
colnames(design) <- c("sample_sizes", "confounder",
                      "Fastsurfer-RCOT", "Fastsurfer-WALD", "Fastsurfer-CPT_KPC", "Fastsurfer-FCIT", "Fastsurfer-CMIknn", "Fastsurfer-PCM",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "Freesurfer-CPT_KPC", "Freesurfer-FCIT", "Freesurfer-CMIknn", "Freesurfer-PCM",
                      "cVAE-RCOT", "cVAE-WALD", "cVAE-CPT_KPC", "cVAE-FCIT", "cVAE-CMIknn", "cVAE-PCM",
                      "MedicalNet-RCOT", "MedicalNet-WALD", "MedicalNet-CPT_KPC", "MedicalNet-FCIT", "MedicalNet-CMIknn", "MedicalNet-PCM")
custom_order <- c("sample_sizes", "confounder",
                  "Fastsurfer-RCOT", "Freesurfer-RCOT","cVAE-RCOT", "MedicalNet-RCOT",
                  "Fastsurfer-CPT_KPC", "Freesurfer-CPT_KPC", "cVAE-CPT_KPC", "MedicalNet-CPT_KPC",
                  "Fastsurfer-FCIT", "Freesurfer-FCIT", "cVAE-FCIT", "MedicalNet-FCIT",
                  "Fastsurfer-CMIknn","Freesurfer-CMIknn", "cVAE-CMIknn", "MedicalNet-CMIknn",
                  "Fastsurfer-PCM","Freesurfer-PCM", "cVAE-PCM", "MedicalNet-PCM",
                  "Fastsurfer-WALD", "Freesurfer-WALD","cVAE-WALD", "MedicalNet-WALD")
#resort columns
design <- design[, custom_order]
design_conf_dim_ci <- design
#design_conf_dim_no_ci <- design



### 2) Data preparation for nested loop over
# loop 1: confounder functional relationship g_z ("linear", "squared", "realistic")
# loop 2: sample size
# fixed confounder dimension 6
# Deep-WALD and Deep-RCOT for all four embedding maps (freesurfer, fastsurfer, condVAE, noisy fastsurfer)
cit_patterns <- "WALD|RCOT|kpc_graph|FCIT|CMIknn|comets_pcm"
cit_files <- grep(cit_patterns, all_files, value=TRUE)
files_conf_relation <- grep("ukb_z6", cit_files, value=TRUE)

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
               confounder = c('linear', 'squared', 'realistic'))#, 'breakpoint3'))
design <- expand.grid(params)
design <- design %>%
  mutate("fastsurfer_fastsurfer-RCOT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-WALD" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-kpc_graph" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-FCIT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-CMIknn" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-comets_pcm" = rep(NA, nrow(design)),
         "freesurfer-RCOT" = rep(NA, nrow(design)),
         "freesurfer-WALD" = rep(NA, nrow(design)),
         "freesurfer-kpc_graph" = rep(NA, nrow(design)),
         "freesurfer-FCIT" = rep(NA, nrow(design)),
         "freesurfer-CMIknn" = rep(NA, nrow(design)),
         "freesurfer-comets_pcm" = rep(NA, nrow(design)),
         "condVAE-RCOT" = rep(NA, nrow(design)),
         "condVAE-WALD" = rep(NA, nrow(design)),
         "condVAE-kpc_graph" = rep(NA, nrow(design)),
         "condVAE-FCIT" = rep(NA, nrow(design)),
         "condVAE-CMIknn" = rep(NA, nrow(design)),
         "condVAE-comets_pcm" = rep(NA, nrow(design)),
         "medicalnet-RCOT" = rep(NA, nrow(design)),
         "medicalnet-WALD" = rep(NA, nrow(design)),
         "medicalnet-kpc_graph" = rep(NA, nrow(design)),
         "medicalnet-FCIT" = rep(NA, nrow(design)),
         "medicalnet-CMIknn" = rep(NA, nrow(design)),
         "medicalnet-comets_pcm" = rep(NA, nrow(design)))

embedding_maps <- c('fastsurfer_fastsurfer', 'freesurfer', 'condVAE', 'medicalnet')
dncits <- c('RCOT', 'WALD', 'kpc_graph', 'FCIT', 'CMIknn', 'comets_pcm')
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

design$confounder <- rep(c('linear', 'squared', 'complex'), each=10)
colnames(design) <- c("sample_sizes", "confounder",
                      "Fastsurfer-RCOT", "Fastsurfer-WALD", "Fastsurfer-CPT_KPC", "Fastsurfer-FCIT", "Fastsurfer-CMIknn", "Fastsurfer-PCM",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "Freesurfer-CPT_KPC", "Freesurfer-FCIT", "Freesurfer-CMIknn", "Freesurfer-PCM",
                      "cVAE-RCOT", "cVAE-WALD", "cVAE-CPT_KPC", "cVAE-FCIT", "cVAE-CMIknn", "cVAE-PCM",
                      "MedicalNet-RCOT", "MedicalNet-WALD", "MedicalNet-CPT_KPC", "MedicalNet-FCIT", "MedicalNet-CMIknn", "MedicalNet-PCM")
custom_order <- c("sample_sizes", "confounder",
                  "Fastsurfer-RCOT", "Freesurfer-RCOT","cVAE-RCOT", "MedicalNet-RCOT",
                  "Fastsurfer-CPT_KPC", "Freesurfer-CPT_KPC", "cVAE-CPT_KPC", "MedicalNet-CPT_KPC",
                  "Fastsurfer-FCIT", "Freesurfer-FCIT", "cVAE-FCIT", "MedicalNet-FCIT",
                  "Fastsurfer-CMIknn","Freesurfer-CMIknn", "cVAE-CMIknn", "MedicalNet-CMIknn",
                  "Fastsurfer-PCM","Freesurfer-PCM", "cVAE-PCM", "MedicalNet-PCM",
                  "Fastsurfer-WALD", "Freesurfer-WALD","cVAE-WALD", "MedicalNet-WALD")
#resort columns
design <- design[, custom_order]
design_conf_relation_ci <- design
#design_conf_relation_no_ci <- design



##### nested loop plots
## Power and T1E together
## conf relation
design <- rbind(design_conf_relation_ci, design_conf_relation_no_ci)
design$'Setting' <- rep(c('No', 'Yes'), each=30)
design$confounder <- rep(rep(c(1,2,3),each=10),2)
## APPENDIX
# without PCM
#design <- design %>% select(-contains("PCM"))
p_conf_relation <- looplot::nested_loop_plot(resdf = design,
                                             x = "sample_sizes",
                                             grid_rows = 'Setting',
                                             steps = "confounder",
                                             #methods = methods_depicted,
                                             steps_y_base = -0.1, steps_y_height = 0.05,,
                                             #legend_breaks = methods_depicted,  # Specify the desired order
                                             #legend_labels = methods_depicted,  # Custom labels if needed
                                             x_name = "Sample size", y_name = "Rejection rate",
                                             spu_x_shift = 1,
                                             colors = palet_discrete[rep(c(1,2,3,4,7), each=3)],
                                             line_linetypes = c(6,1,3),
                                             point_size = 3,
                                             line_size = 1.5,
                                             point_shapes = c(17,19,15),
                                             steps_values_annotate = TRUE, steps_annotation_size = 6, steps_color='grey31',steps_annotation_color='grey31',
                                             hline_intercept = c(0,0.05),
                                             y_expand_add = c(0.1,0.15),
                                             line_alpha =0.6,
                                             point_alpha = 0.8,
                                             legend_name = "DNCIT",
                                             base_size = 24,
                                             replace_labels = list(
                                               Setting = c('No'='T1E',
                                                           'Yes'='Power'),
                                               confounder = c('1'='linear',
                                                              '2'='squared',
                                                              '3'='complex')
                                             ),
                                             grid_labeller = labeller('No'='T1E', 'Yes'='Power'),
                                             post_processing = list(
                                               add_custom_theme = list(
                                                 axis.text.x = ggplot2::element_text(angle = -90,
                                                                                     vjust = 0.5,
                                                                                     size = 15)
                                               )
                                             ))
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_rel_power_T1E_app.png'), p_conf_relation, width = 16, height = 16, dpi = 300)
print(p_conf_relation)

## MAIN TEXT (3 columns)
design_main_text <- design %>%
  mutate(across(contains("Fastsurfer-") ,#| contains("cVAE-"),
                ~ ifelse(Setting == "No", NA, .)))
#remove tests without power from power plots
set_to_NA_because_no_power <- c()
for(col in colnames(design_main_text)[3:(ncol(design_main_text)-1)]){
  for(confounder in c(1,2,3)){
    if(grepl("CMIknn", col)){
      max_sample_size <- 1100
    }else{
      max_sample_size <- 10000
    }
    if(any(design_main_text$confounder == confounder & design_main_text$sample_sizes == max_sample_size & design_main_text[[col]] < 0.15 & design_main_text$Setting == "Yes", na.rm=TRUE)){
      print(col)
      print(confounder)
      design_main_text[which(design_main_text$confounder == confounder & design_main_text$Setting == "Yes"), col] <- NA
      set_to_NA_because_no_power <- c(set_to_NA_because_no_power, paste(col, confounder))
    }
    for (sample_size in c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)){
      if(any(design_main_text$confounder == confounder & design_main_text$sample_sizes == sample_size & design_main_text[[col]] > 0.15 & design_main_text$Setting == "No", na.rm=TRUE)){
        print(col)
        print(confounder)
        design_main_text[which(design_main_text$confounder == confounder & design_main_text$sample_sizes == sample_size & design_main_text$Setting == "Yes"), col] <- NA
        set_to_NA_because_no_power <- c(set_to_NA_because_no_power, paste(col, confounder))
      }
    }
  }
}

# without PCM
# design_main_text <- design_main_text %>% select(-contains("PCM"))
#withour Fastsurfer
design_main_text <- design_main_text %>% select(-contains("Fastsurfer"))
methods_depicted <- colnames(design_main_text)[-c(1:2, ncol(design_main_text))]

p_conf_relation <- looplot::nested_loop_plot(resdf = design_main_text,
                                             x = "sample_sizes",
                                             grid_rows = 'Setting',
                                             steps = "confounder",
                                             methods = methods_depicted,
                                             steps_y_base = -0.1, steps_y_height = 0.05,,
                                             legend_breaks = methods_depicted,  # Specify the desired order
                                             legend_labels = methods_depicted,  # Custom labels if needed
                                             x_name = "Sample size", y_name = "Rejection rate",
                                             spu_x_shift = 1,
                                             colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                             line_linetypes = c(1,2,3),
                                             point_size = 3,
                                             line_size = 1.5,
                                             point_shapes = c(19,17,15),
                                             steps_values_annotate = TRUE, steps_annotation_size = 6,
                                             hline_intercept = c(0,0.05),
                                             y_expand_add = c(0.1,0.15),
                                             line_alpha =0.6,
                                             point_alpha = 0.8,
                                             legend_name = "DNCIT",
                                             base_size = 24,
                                             replace_labels = list(
                                               Setting = c('No'='T1E',
                                                           'Yes'='Power'),
                                               confounder = c('1'='linear',
                                                              '2'='squared',
                                                              '3'='complex')
                                             ),
                                             grid_labeller = labeller('No'='T1E', 'Yes'='Power'),
                                             post_processing = list(
                                               add_custom_theme = list(
                                                 axis.text.x = ggplot2::element_text(angle = -90,
                                                                                     vjust = 0.5,
                                                                                     size = 15)
                                               )
                                             ))
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_rel_power_T1E.png'), p_conf_relation, width = 16, height = 10, dpi = 300)
print(p_conf_relation)

## Power and T1E separate
## Appendix
# CI
design_ci <- design[!(design$Setting=='Yes'),]
p_conf_relation_ci <- looplot::nested_loop_plot(resdf = design_ci,
                                                x = "sample_sizes",
                                                grid_rows = 'Setting',
                                                steps = "confounder",
                                                methods = methods_depicted,
                                                steps_y_base = -0.15/2, steps_y_height = 0.05/2,
                                                legend_breaks = methods_depicted,  # Specify the desired order
                                                legend_labels = methods_depicted,  # Custom labels if needed
                                                x_name = "Sample size", y_name = "Rejection rate",
                                                spu_x_shift = 1,
                                                colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                                line_linetypes = c(1,2,3),
                                                point_size = 4,
                                                line_size = 1.5,
                                                point_shapes = c(19,17,15),
                                                steps_values_annotate = TRUE, steps_annotation_size = 6, steps_color='grey31',steps_annotation_color='grey31',
                                                hline_intercept = c(0,0.05),
                                                hline_linetype =c(1),
                                                hline_size = c(0.5,1.5),
                                                hline_colour = "black",
                                                y_expand_add = c(0.1/2,0.15/2),
                                                line_alpha =0.6,
                                                point_alpha = 0.8,
                                                legend_name = "DNCIT",
                                                base_size = 24,
                                                replace_labels = list(
                                                  Setting = c('No'='T1E'),
                                                  confounder = c('1'='linear',
                                                                 '2'='squared',
                                                                 '3'='complex')
                                                ),
                                                grid_labeller = labeller('No'='T1E'),
                                                post_processing = list(
                                                  add_custom_theme = list(
                                                    axis.text.x = ggplot2::element_text(angle = -90,
                                                                                        vjust = 0.5,
                                                                                        size = 15)
                                                  )
                                                ))
##No CI
design_no_ci <- design[!(design$Setting=='No'),]
p_conf_relation_no_ci <- looplot::nested_loop_plot(resdf = design_no_ci,
                                                   x = "sample_sizes",
                                                   grid_rows = 'Setting',
                                                   steps = "confounder",
                                                   methods = methods_depicted,
                                                   steps_y_base = -0.15, steps_y_height = 0.05,
                                                   legend_breaks = methods_depicted,  # Specify the desired order
                                                   legend_labels = methods_depicted,  # Custom labels if needed
                                                   x_name = "Sample size", y_name = "Rejection rate",
                                                   spu_x_shift = 1,
                                                   colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                                   line_linetypes = c(1,2,3),
                                                   point_size = 4,
                                                   line_size = 1.5,
                                                   point_shapes = c(19,17,15),
                                                   steps_values_annotate = TRUE, steps_annotation_size = 6, steps_color='grey31',steps_annotation_color='grey31',
                                                   hline_intercept = c(0),
                                                   hline_linetype =1,
                                                   hline_size = c(0.5),
                                                   hline_colour = "black",
                                                   y_expand_add = c(0.1,0.15),
                                                   y_breaks = seq(0,1,0.2),
                                                   line_alpha =0.6,
                                                   point_alpha = 0.8,
                                                   legend_name = "DNCIT",
                                                   base_size = 24,
                                                   replace_labels = list(
                                                     Setting = c(#'No'='T1E',
                                                       'Yes'='Power'),
                                                     confounder = c('1'='linear',
                                                                    '2'='squared',
                                                                    '3'='complex')
                                                   ),
                                                   grid_labeller = labeller(#'No'='T1E',
                                                     'Yes'='Power'),
                                                   post_processing = list(
                                                     add_custom_theme = list(
                                                       axis.text.x = ggplot2::element_text(angle = -90,
                                                                                           vjust = 0.5,
                                                                                           size = 15)
                                                     )
                                                   ))

# legend
design_legend <- design %>%
  mutate(across(contains("Fastsurfer-"),# | contains("cVAE-"),
                ~ ifelse(Setting == "No", NA, .)))
# without PCM
#design_legend <- design_legend %>% select(-contains("PCM"))
#withour Fastsurfer
design_legend <- design_legend %>% select(-contains("Fastsurfer"))
methods_depicted <- colnames(design_legend)[-c(1:2, ncol(design_legend))]
p_conf_relation_legend <- looplot::nested_loop_plot(resdf = design_legend,
                                                    x = "sample_sizes",
                                                    grid_rows = 'Setting',
                                                    steps = "confounder",
                                                    methods = methods_depicted,
                                                    steps_y_base = -0.1, steps_y_height = 0.05,,
                                                    legend_breaks = methods_depicted,  # Specify the desired order
                                                    legend_labels = methods_depicted,  # Custom labels if needed
                                                    x_name = "Sample size", y_name = "Rejection rate",
                                                    spu_x_shift = 1,
                                                    colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                                    line_linetypes = c(1,2,3),
                                                    point_size = 4,
                                                    line_size = 1.5,
                                                    point_shapes = c(19,17,15),
                                                    steps_values_annotate = TRUE, steps_annotation_size = 6, steps_color='grey31',steps_annotation_color='grey31',
                                                    hline_intercept = c(0,0.05),
                                                    y_expand_add = c(0.1,0.15/2),
                                                    y_breaks = seq(0,1,0.2),
                                                    line_alpha =0.6,
                                                    point_alpha = 0.8,
                                                    legend_name = "DNCIT",
                                                    base_size = 24,
                                                    replace_labels = list(
                                                      Setting = c(#'No'='T1E',
                                                        'Yes'='Power'),
                                                      confounder = c('1'='linear',
                                                                     '2'='squared',
                                                                     '3'='complex')
                                                    ),
                                                    grid_labeller = labeller(#'No'='T1E',
                                                      'Yes'='Power'),
                                                    post_processing = list(
                                                      add_custom_theme = list(
                                                        axis.text.x = ggplot2::element_text(angle = -90,
                                                                                            vjust = 0.5,
                                                                                            size = 15)
                                                      )
                                                    ))
p_conf_relation_ci_mod <- p_conf_relation_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),  # Remove x-axis label from first plot
        axis.text.x = element_blank(),   # Remove x-axis ticks from first plot
        axis.ticks.x = element_blank())
p_conf_relation_no_ci_mod <- p_conf_relation_no_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "cm"),  # Increase overall size of the legend symbols
        legend.spacing.x = unit(0.5, "cm"), # Increase horizontal spacing between legend items
        legend.spacing.y = unit(0.5, "cm")) + # Adjust vertical spacing if needed
  guides(
    colour = guide_legend(nrow = 2, keywidth = 2, keyheight = 2),  # Increase width and height of the legend symbols
    shape = guide_legend(nrow = 2, keywidth = 2, keyheight = 2)    # Same for shape legend if needed
  )

combined_plot <- plot_grid(
  p_conf_relation_ci_mod,  # First plot
  p_conf_relation_no_ci_mod,  # Second plot with x-axis labels
  ncol = 1,  # Stack them vertically
  align = "v",  # Align the vertical axes
  rel_heights = c(1, 1.2)  # Adjust relative heights if needed
)
y_axis_label <- ggdraw() +
  draw_label("Rejection Rate", x = 0.5, y = 0.6, angle = 90, vjust = 0.5, size = 24) +
  theme(panel.background = element_rect(fill = "white", colour = NA), 
        plot.background = element_rect(fill = "white", colour = NA))
combined_plot <- plot_grid(
  y_axis_label,  # Y-axis label on the left
  combined_plot,  # The stacked plots
  ncol = 2,  # Arrange in two columns (label + plots)
  rel_widths = c(0.05, 1)  # Give space for the y-axis label
)
# Apply guides() to modify the legend appearance (e.g., number of rows, key size)
p_conf_relation_no_ci_legend <- p_conf_relation_legend +  # Customize the legend
  guides(colour = guide_legend(nrow = 3, keywidth = 1.5, keyheight=1.2, override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 17),
        legend.title = element_text(size = 18))

# Use cowplot's get_legend() to extract the modified legend
legend <- get_legend(p_conf_relation_no_ci_legend)
p_conf_relation <- plot_grid(
  combined_plot,  # Stacked plots
  legend,  # Legend below the plots
  ncol = 1,  # Legend below, so keep 1 column
  rel_heights = c(1, 0.1)  # Adjust height ratios if needed
)
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_rel_power_T1E_app.png'), p_conf_relation, width = 16, height = 16, dpi = 300)
print(p_conf_relation)


## MAIN TEXT (2 columns)
# CI
design_main_text_ci <- design_main_text[!(design_main_text$confounder ==2 | design_main_text$Setting=='Yes'),] %>% 
    select(-contains("cVAE-") & -contains("MedicalNet-"))
methods_depicted_ci <- colnames(design_main_text_ci)[-c(1:2, ncol(design_main_text_ci))]
p_conf_relation_ci <- looplot::nested_loop_plot(resdf = design_main_text_ci,
                                             x = "sample_sizes",
                                             grid_rows = 'Setting',
                                             steps = "confounder",
                                             methods = methods_depicted_ci,
                                             steps_y_base = -0.15/2, steps_y_height = 0.05/2,
                                             legend_breaks = methods_depicted_ci,  # Specify the desired order
                                             legend_labels = methods_depicted_ci,  # Custom labels if needed
                                             x_name = "Sample size", y_name = "Rejection rate",
                                             spu_x_shift = 1,
                                             colors = palet_discrete[rep(c(1,2,7,4,3,6))],#, each=3)],
                                             line_linetypes = c(1),#,2,3),
                                             point_size = 4,
                                             line_size = 1.5,
                                             point_shapes = c(19),#,17,15),
                                             steps_values_annotate = TRUE, steps_annotation_size = 8, steps_color='grey31',steps_annotation_color='grey31',
                                             hline_intercept = c(0,0.05),
                                             hline_linetype =c(1),
                                             hline_size = c(0.5,1.5),
                                             hline_colour = "black",
                                             y_expand_add = c(0.15,0.15),
                                             line_alpha =0.6,
                                             point_alpha = 0.8,
                                             legend_name = "DNCIT",
                                             base_size = 30,
                                             replace_labels = list(
                                               Setting = c('No'='T1E'),
                                               confounder = c('1'='linear',
                                                              #'2'='squared',
                                                              '3'='complex')
                                             ),
                                             grid_labeller = labeller('No'='T1E'),
                                             post_processing = list(
                                               add_custom_theme = list(
                                                 axis.text.x = ggplot2::element_text(angle = -90,
                                                                                     vjust = 0.5)
                                               )
                                             ))
##No CI
design_main_text_no_ci <- design_main_text[!(design_main_text$confounder ==2 | design_main_text$Setting=='No'),]
p_conf_relation_no_ci <- looplot::nested_loop_plot(resdf = design_main_text_no_ci,
                                             x = "sample_sizes",
                                             grid_rows = 'Setting',
                                             steps = "confounder",
                                             methods = methods_depicted,
                                             steps_y_base = -0.15, steps_y_height = 0.05,
                                             legend_breaks = methods_depicted,  # Specify the desired order
                                             legend_labels = methods_depicted,  # Custom labels if needed
                                             x_name = "Sample size", y_name = "Rejection rate",
                                             spu_x_shift = 1,
                                             colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                             line_linetypes = c(1,2,3),
                                             point_size = 4,
                                             line_size = 1.5,
                                             point_shapes = c(19,17,15),
                                             steps_values_annotate = TRUE, steps_annotation_size = 8, steps_color='grey31',steps_annotation_color='grey31',
                                             hline_intercept = c(0),
                                             hline_linetype =1,
                                             hline_size = c(0.5),
                                             hline_colour = "black",
                                             y_expand_add = c(0.1,0.15/2),
                                             y_breaks = seq(0,1,0.2),
                                             line_alpha =0.6,
                                             point_alpha = 0.8,
                                             legend_name = "DNCIT",
                                             base_size = 30,
                                             replace_labels = list(
                                               Setting = c(#'No'='T1E',
                                                 'Yes'='Power'),
                                               confounder = c('1'='linear',
                                                              #'2'='squared',
                                                              '3'='complex')
                                             ),
                                             grid_labeller = labeller(#'No'='T1E',
                                               'Yes'='Power'),
                                             post_processing = list(
                                               add_custom_theme = list(
                                                 axis.text.x = ggplot2::element_text(angle = -90,
                                                                                     vjust = 0.5)
                                               )
                                             ))

# legend
design_legend <- design %>%
  mutate(across(contains("Fastsurfer-"), # | contains("cVAE-"),
                ~ ifelse(Setting == "No", NA, .)))
# without PCM
# design_legend <- design_legend %>% select(-contains("PCM"))
#withour Fastsurfer
design_legend <- design_legend %>% select(-contains("Fastsurfer"))
methods_depicted <- colnames(design_legend)[-c(1:2, ncol(design_legend))]
p_conf_relation_legend <- looplot::nested_loop_plot(resdf = design_legend,
                                                   x = "sample_sizes",
                                                   grid_rows = 'Setting',
                                                   steps = "confounder",
                                                   methods = methods_depicted,
                                                   steps_y_base = -0.1, steps_y_height = 0.05,,
                                                   legend_breaks = methods_depicted,  # Specify the desired order
                                                   legend_labels = methods_depicted,  # Custom labels if needed
                                                   x_name = "Sample size", y_name = "Rejection rate",
                                                   spu_x_shift = 1,
                                                   colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                                   line_linetypes = c(1,2,3),
                                                   point_size = 4,
                                                   line_size = 1.5,
                                                   point_shapes = c(19,17,15),
                                                   steps_values_annotate = TRUE, steps_annotation_size = 7, steps_color='grey31',steps_annotation_color='grey31',
                                                   hline_intercept = c(0,0.05),
                                                   y_expand_add = c(0.1,0.15/2),
                                                   y_breaks = seq(0,1,0.2),
                                                   line_alpha =0.6,
                                                   point_alpha = 0.8,
                                                   legend_name = "DNCIT",
                                                   base_size = 28,
                                                   replace_labels = list(
                                                     Setting = c(#'No'='T1E',
                                                       'Yes'='Power'),
                                                     confounder = c('1'='linear',
                                                                    #'2'='squared',
                                                                    '3'='complex')
                                                   ),
                                                   grid_labeller = labeller(#'No'='T1E',
                                                     'Yes'='Power'),
                                                   post_processing = list(
                                                     add_custom_theme = list(
                                                       axis.text.x = ggplot2::element_text(angle = -90,
                                                                                           vjust = 0.5)
                                                     )
                                                   ))
p_conf_relation_ci_mod <- p_conf_relation_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),  # Remove x-axis label from first plot
        axis.text.x = element_blank(),   # Remove x-axis ticks from first plot
        axis.ticks.x = element_blank())
p_conf_relation_no_ci_mod <- p_conf_relation_no_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "cm"),  # Increase overall size of the legend symbols
        legend.spacing.x = unit(0.5, "cm"), # Increase horizontal spacing between legend items
        legend.spacing.y = unit(0.5, "cm")) + # Adjust vertical spacing if needed
  guides(
    colour = guide_legend(nrow = 2, keywidth = 2, keyheight = 2),  # Increase width and height of the legend symbols
    shape = guide_legend(nrow = 2, keywidth = 2, keyheight = 2)    # Same for shape legend if needed
  )

combined_plot <- plot_grid(
  p_conf_relation_ci_mod,  # First plot
  p_conf_relation_no_ci_mod,  # Second plot with x-axis labels
  ncol = 1,  # Stack them vertically
  align = "v",  # Align the vertical axes
  rel_heights = c(1, 1.2)  # Adjust relative heights if needed
)
y_axis_label <- ggdraw() +
  draw_label("Rejection Rate", x = 0.5, y = 0.6, angle = 90, vjust = 0.5, size = 24) +
  theme(panel.background = element_rect(fill = "white", colour = NA), 
        plot.background = element_rect(fill = "white", colour = NA))
combined_plot <- plot_grid(
  y_axis_label,  # Y-axis label on the left
  combined_plot,  # The stacked plots
  ncol = 2,  # Arrange in two columns (label + plots)
  rel_widths = c(0.05, 1)  # Give space for the y-axis label
)
# Apply guides() to modify the legend appearance (e.g., number of rows, key size)
p_conf_relation_no_ci_legend <- p_conf_relation_legend +  # Customize the legend
  guides(colour = guide_legend(nrow = 3, keywidth = 1.5, keyheight=1.2, override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 19.3),
        legend.title = element_text(size = 20))

# Use cowplot's get_legend() to extract the modified legend
legend <- get_legend(p_conf_relation_no_ci_legend)
p_conf_relation <- plot_grid(
  combined_plot,  # Stacked plots
  legend,  # Legend below the plots
  ncol = 1,  # Legend below, so keep 1 column
  rel_heights = c(1, 0.1)  # Adjust height ratios if needed
)
# ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_rel_power_T1E.png'), p_conf_relation, width = 18, height = 16, dpi = 300)
print(p_conf_relation)

############ conf dimension
design <- rbind(design_conf_dim_ci, design_conf_dim_no_ci)
design$'Setting' <- rep(c('No', 'Yes'), each=60)
# confounder dimension (1,2,4,6,10,15) [APPENDIX] vs confounder dimension (1,2,10) [MAIN TEXT]
# APPENDIX
# without PCM
#design <- design %>% select(-contains("PCM"))
p_conf_dim = looplot::nested_loop_plot(resdf = design,
                                       x = "sample_sizes",
                                       grid_rows = 'Setting',
                                       steps  = "confounder",
                                       steps_y_base = -0.1, steps_y_height = 0.05,
                                       x_name = "Sample size", y_name = "Rejection rate",
                                       spu_x_shift = 1,
                                       colors = palet_discrete[c(rep(1,3), rep(2,3), rep(7,3), rep(4,3), rep(3,3), rep(6,3))],
                                       line_linetypes = c(1,2,3),
                                       point_size = 3,
                                       line_size = 1,
                                       point_shapes = c(19,17,15),
                                       steps_values_annotate = TRUE, steps_annotation_size = 6,
                                       y_expand_add = c(0.1,0.15),
                                       hline_intercept = c(0,0.05),
                                       line_alpha =0.6,
                                       point_alpha = 0.8,
                                       legend_name = "DNCIT",
                                       base_size = 22,
                                       replace_labels = list(
                                         Setting = c('No'='T1E',
                                                     'Yes'='Power')
                                       ),
                                       grid_labeller = labeller('No'='T1E', 'Yes'='Power'),
                                       post_processing = list(
                                         add_custom_theme = list(
                                           axis.text.x = ggplot2::element_text(angle = -90,
                                                                               vjust = 0.5,
                                                                               size = 15)
                                         )
                                       ))

#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_dim.png'), p_conf_dim, width = 16, height = 10, dpi = 300)
print(p_conf_dim)

# MAIN TEXT:
# design <- design[design$'confounder' %in% c(1,2,10),]
#set for T1E Fastsurfer an cVAE to NA
design_main_text <- design %>%
  mutate(across(contains("Fastsurfer-") ,#| contains("cVAE-"),
                ~ ifelse(Setting == "No", NA, .)))
#remove tests without power from power plots
set_to_NA_because_no_power <- c()
for(col in colnames(design_main_text)[3:(ncol(design_main_text)-1)]){
  for(confounder in c(1,2,4,6,10,15)){
    if(grepl("CMIknn", col)){
      max_sample_size <- 1100
    }else{
      max_sample_size <- 10000
    }
    if(any(design_main_text$confounder == confounder & design_main_text$sample_sizes == max_sample_size & design_main_text[[col]] < 0.15 & design_main_text$Setting == "Yes", na.rm=TRUE)){
      print(col)
      print(confounder)
      design_main_text[which(design_main_text$confounder == confounder & design_main_text$Setting == "Yes"), col] <- NA
      set_to_NA_because_no_power <- c(set_to_NA_because_no_power, paste(col, confounder))
    }
    for (sample_size in c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)){
      if(any(design_main_text$confounder == confounder & design_main_text$sample_sizes == sample_size & design_main_text[[col]] > 0.15 & design_main_text$Setting == "No", na.rm=TRUE)){
        print(col)
        print(confounder)
        design_main_text[which(design_main_text$confounder == confounder & design_main_text$sample_sizes == sample_size & design_main_text$Setting == "Yes"), col] <- NA
        set_to_NA_because_no_power <- c(set_to_NA_because_no_power, paste(col, confounder))
      }
    }
  }
}

# without PCM
# design_main_text <- design_main_text %>% select(-contains("PCM"))
#withour Fastsurfer
design_main_text <- design_main_text %>% select(-contains("Fastsurfer-"))
design_main_text <- design_main_text %>%
  mutate(across(contains("Fastsurfer-") | contains("cVAE-") | contains("MedicalNet-"),
                ~ ifelse(Setting == "No", NA, .)))
methods_depicted_ci <- colnames(design_main_text)[-c(1:2, ncol(design_main_text))]
p_conf_dim = looplot::nested_loop_plot(resdf = design_main_text,
                                       x = "sample_sizes",
                                       grid_rows = 'Setting',
                                       steps  = "confounder",
                                       methods = methods_depicted_ci,
                                       steps_y_base = -0.1, steps_y_height = 0.05,,
                                       legend_breaks = methods_depicted_ci,  # Specify the desired order
                                       legend_labels = methods_depicted_ci,  # Custom labels if needed
                                       x_name = "Sample size", y_name = "Rejection rate",
                                       spu_x_shift = 1,
                                       colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                       line_linetypes = rep(c(1,2,3),4),
                                       point_size = 3,
                                       line_size = 1,
                                       point_shapes = c(19,17,15),
                                       steps_values_annotate = TRUE, steps_annotation_size = 6,
                                       y_expand_add = c(0.1,0.15),
                                       hline_intercept = c(0,0.05),
                                       line_alpha =0.6,
                                       point_alpha = 0.8,
                                       legend_name = "DNCIT",
                                       base_size = 24,
                                       replace_labels = list(
                                         Setting = c('No'='T1E',
                                                     'Yes'='Power')
                                       ),
                                       grid_labeller = labeller('No'='T1E', 'Yes'='Power'),
                                       post_processing = list(
                                         add_custom_theme = list(
                                           axis.text.x = ggplot2::element_text(angle = -90,
                                                                               vjust = 0.5,
                                                                               size = 15)
                                         )
                                       ))
ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_dim_main_text.png'), p_conf_dim, width = 16, height = 10, dpi = 300)
print(p_conf_dim)

### Power and T1E separate
## APPENDIX
# CI
colnames(design)[2] <- "confounder dimension"
design_ci <- design[!(design$Setting=='Yes'),]
p_conf_dim_ci <- looplot::nested_loop_plot(resdf = design_ci,
                                           x = "sample_sizes",
                                           grid_rows = 'Setting',
                                           steps = "confounder dimension",
                                           methods = methods_depicted,
                                           steps_y_base = -0.15*4/5, steps_y_height = 0.05*4/5,
                                           legend_breaks = methods_depicted,  # Specify the desired order
                                           legend_labels = methods_depicted,  # Custom labels if needed
                                           x_name = "Sample size", y_name = "Rejection rate",
                                           spu_x_shift = 1,
                                           colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                           line_linetypes = c(1,2,3),
                                           point_size = 4,
                                           line_size = 1.5,
                                           point_shapes = c(19,17,15),
                                           steps_values_annotate = TRUE, steps_annotation_size = 6, steps_color='grey31',steps_annotation_color='grey31',
                                           hline_intercept = c(0,0.05),
                                           hline_linetype =c(1),
                                           hline_size = c(0.5,1.5),
                                           hline_colour = "black",
                                           y_expand_add = c(0.2,0.15),
                                           line_alpha =0.6,
                                           point_alpha = 0.8,
                                           legend_name = "DNCIT",
                                           base_size = 24,
                                           replace_labels = list(Setting = c('No'='T1E')
                                           ),
                                           grid_labeller = labeller('No'='T1E'),
                                           post_processing = list(
                                             add_custom_theme = list(
                                               axis.text.x = ggplot2::element_text(angle = -90,
                                                                                   vjust = 0.5,
                                                                                   size = 15)
                                             )
                                           ))
##No CI
design_no_ci <- design[!(design$Setting=='No'),]
p_conf_dim_no_ci <- looplot::nested_loop_plot(resdf = design_no_ci,
                                              x = "sample_sizes",
                                              grid_rows = 'Setting',
                                              steps = "confounder dimension",
                                              methods = methods_depicted,
                                              steps_y_base = -0.15, steps_y_height = 0.05,
                                              legend_breaks = methods_depicted,  # Specify the desired order
                                              legend_labels = methods_depicted,  # Custom labels if needed
                                              x_name = "Sample size", y_name = "Rejection rate",
                                              spu_x_shift = 1,
                                              colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                              line_linetypes = c(1,2,3),
                                              point_size = 4,
                                              line_size = 1.5,
                                              point_shapes = c(19,17,15),
                                              steps_values_annotate = TRUE, steps_annotation_size = 6, steps_color='grey31',steps_annotation_color='grey31',
                                              hline_intercept = c(0),
                                              hline_linetype =1,
                                              hline_size = c(0.5),
                                              hline_colour = "black",
                                              y_expand_add = c(0.2,0.15),
                                              y_breaks = seq(0,1,0.2),
                                              line_alpha =0.6,
                                              point_alpha = 0.8,
                                              legend_name = "DNCIT",
                                              base_size = 24,
                                              replace_labels = list(
                                                Setting = c('Yes'='Power')
                                              ),
                                              grid_labeller = labeller('Yes'='Power'),
                                              post_processing = list(
                                                add_custom_theme = list(
                                                  axis.text.x = ggplot2::element_text(angle = -90,
                                                                                      vjust = 0.5,
                                                                                      size = 15)
                                                )
                                              ))

# legend
design_legend <- design %>%
  mutate(across(contains("Fastsurfer-"),# | contains("cVAE-"),
                ~ ifelse(Setting == "No", NA, .)))
# without PCM
#design_legend <- design_legend %>% select(-contains("PCM"))
#withour Fastsurfer
design_legend <- design_legend %>% select(-contains("Fastsurfer"))
methods_depicted <- colnames(design_legend)[-c(1:2, ncol(design_legend))]
p_conf_dim_legend <- looplot::nested_loop_plot(resdf = design_legend,
                                               x = "sample_sizes",
                                               grid_rows = 'Setting',
                                               steps = "confounder dimension",
                                               methods = methods_depicted,
                                               steps_y_base = -0.1, steps_y_height = 0.05,,
                                               legend_breaks = methods_depicted,  # Specify the desired order
                                               legend_labels = methods_depicted,  # Custom labels if needed
                                               x_name = "Sample size", y_name = "Rejection rate",
                                               spu_x_shift = 1,
                                               colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                               line_linetypes = c(1,2,3),
                                               point_size = 4,
                                               line_size = 1.5,
                                               point_shapes = c(19,17,15),
                                               steps_values_annotate = TRUE, steps_annotation_size = 6, steps_color='grey31',steps_annotation_color='grey31',
                                               hline_intercept = c(0,0.05),
                                               y_expand_add = c(0.1,0.15/2),
                                               y_breaks = seq(0,1,0.2),
                                               line_alpha =0.6,
                                               point_alpha = 0.8,
                                               legend_name = "DNCIT",
                                               base_size = 24,
                                               replace_labels = list(
                                                 Setting = c(#'No'='T1E',
                                                   'Yes'='Power')
                                               ),
                                               grid_labeller = labeller(#'No'='T1E',
                                                 'Yes'='Power'),
                                               post_processing = list(
                                                 add_custom_theme = list(
                                                   axis.text.x = ggplot2::element_text(angle = -90,
                                                                                       vjust = 0.5,
                                                                                       size = 15)
                                                 )
                                               ))
p_conf_dim_ci_mod <- p_conf_dim_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),  # Remove x-axis label from first plot
        axis.text.x = element_blank(),   # Remove x-axis ticks from first plot
        axis.ticks.x = element_blank())
p_conf_dim_no_ci_mod <- p_conf_dim_no_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "cm"),  # Increase overall size of the legend symbols
        legend.spacing.x = unit(0.5, "cm"), # Increase horizontal spacing between legend items
        legend.spacing.y = unit(0.5, "cm")) + # Adjust vertical spacing if needed
  guides(
    colour = guide_legend(nrow = 2, keywidth = 2, keyheight = 2),  # Increase width and height of the legend symbols
    shape = guide_legend(nrow = 2, keywidth = 2, keyheight = 2)    # Same for shape legend if needed
  )

combined_plot <- plot_grid(
  p_conf_dim_ci_mod,  # First plot
  p_conf_dim_no_ci_mod,  # Second plot with x-axis labels
  ncol = 1,  # Stack them vertically
  align = "v",  # Align the vertical axes
  rel_heights = c(1, 1.2)  # Adjust relative heights if needed
)
y_axis_label <- ggdraw() +
  draw_label("Rejection Rate", x = 0.5, y = 0.6, angle = 90, vjust = 0.5, size = 24) +
  theme(panel.background = element_rect(fill = "white", colour = NA), 
        plot.background = element_rect(fill = "white", colour = NA))
combined_plot <- plot_grid(
  y_axis_label,  # Y-axis label on the left
  combined_plot,  # The stacked plots
  ncol = 2,  # Arrange in two columns (label + plots)
  rel_widths = c(0.05, 1)  # Give space for the y-axis label
)
# Apply guides() to modify the legend appearance (e.g., number of rows, key size)
p_conf_dim_no_ci_legend <- p_conf_dim_legend +  # Customize the legend
  guides(colour = guide_legend(nrow = 3, keywidth = 2, keyheight=1.2, override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(size = 20))
# Use cowplot's get_legend() to extract the modified legend
legend <- get_legend(p_conf_dim_no_ci_legend)
p_conf_dim <- plot_grid(
  combined_plot,  # Stacked plots
  legend,  # Legend below the plots
  ncol = 1,  # Legend below, so keep 1 column
  rel_heights = c(1, 0.15)  # Adjust height ratios if needed
)
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_dim.png'), p_conf_dim, width = 18, height = 15, dpi = 300)
print(p_conf_dim)


## MAIN TEXT (2 columns)
# CI
design_main_text <- design_main_text[design_main_text$'confounder' %in% c(1,10),]
colnames(design_main_text)[2] <- "confounder dimension"
design_main_text_ci <- design_main_text[!(design_main_text$confounder ==2 | design_main_text$Setting=='Yes'),]
p_conf_dim_ci <- looplot::nested_loop_plot(resdf = design_main_text_ci,
                                                x = "sample_sizes",
                                                grid_rows = 'Setting',
                                                steps = "confounder dimension",
                                                methods = methods_depicted,
                                                steps_y_base = -0.15*4/5, steps_y_height = 0.05*4/5,
                                                legend_breaks = methods_depicted,  # Specify the desired order
                                                legend_labels = methods_depicted,  # Custom labels if needed
                                                x_name = "Sample size", y_name = "Rejection rate",
                                                spu_x_shift = 1,
                                                colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                                line_linetypes = c(1,2,3),
                                                point_size = 4,
                                                line_size = 1.5,
                                                point_shapes = c(19,17,15),
                                                steps_values_annotate = TRUE, steps_annotation_size = 8, steps_color='grey31',steps_annotation_color='grey31',
                                                hline_intercept = c(0,0.05),
                                                hline_linetype =c(1),
                                                hline_size = c(0.5,1.5),
                                                hline_colour = "black",
                                                y_expand_add = c(0.1*4/5,0.15*4/5),
                                                line_alpha =0.6,
                                                point_alpha = 0.8,
                                                legend_name = "DNCIT",
                                                base_size = 30,
                                                replace_labels = list(Setting = c('No'='T1E')
                                                ),
                                                grid_labeller = labeller('No'='T1E'),
                                                post_processing = list(
                                                  add_custom_theme = list(
                                                    axis.text.x = ggplot2::element_text(angle = -90,
                                                                                        vjust = 0.5)
                                                  )
                                                ))
##No CI
design_main_text_no_ci <- design_main_text[!(design_main_text$confounder ==2 | design_main_text$Setting=='No'),]
p_conf_dim_no_ci <- looplot::nested_loop_plot(resdf = design_main_text_no_ci,
                                                   x = "sample_sizes",
                                                   grid_rows = 'Setting',
                                                   steps = "confounder dimension",
                                                   methods = methods_depicted,
                                                   steps_y_base = -0.15, steps_y_height = 0.05,
                                                   legend_breaks = methods_depicted,  # Specify the desired order
                                                   legend_labels = methods_depicted,  # Custom labels if needed
                                                   x_name = "Sample size", y_name = "Rejection rate",
                                                   spu_x_shift = 1,
                                                   colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                                   line_linetypes = c(1,2,3),
                                                   point_size = 4,
                                                   line_size = 1.5,
                                                   point_shapes = c(19,17,15),
                                                   steps_values_annotate = TRUE, steps_annotation_size = 8, steps_color='grey31',steps_annotation_color='grey31',
                                                   hline_intercept = c(0),
                                                   hline_linetype =1,
                                                   hline_size = c(0.5),
                                                   hline_colour = "black",
                                                   y_expand_add = c(0.1,0.15),
                                                   y_breaks = seq(0,1,0.2),
                                                   line_alpha =0.6,
                                                   point_alpha = 0.8,
                                                   legend_name = "DNCIT",
                                                   base_size = 30,
                                                   replace_labels = list(
                                                     Setting = c('Yes'='Power')
                                                   ),
                                                   grid_labeller = labeller('Yes'='Power'),
                                                   post_processing = list(
                                                     add_custom_theme = list(
                                                       axis.text.x = ggplot2::element_text(angle = -90,
                                                                                           vjust = 0.5)
                                                     )
                                                   ))

# legend
design_legend <- design %>%
  mutate(across(contains("Fastsurfer-"),# | contains("cVAE-"),
                ~ ifelse(Setting == "No", NA, .)))
# without PCM
design_legend <- design_legend %>% select(-contains("PCM"))
#withour Fastsurfer
design_legend <- design_legend %>% select(-contains("Fastsurfer"))
methods_depicted <- colnames(design_legend)[-c(1:2, ncol(design_legend))]
p_conf_dim_legend <- looplot::nested_loop_plot(resdf = design_legend,
                                                    x = "sample_sizes",
                                                    grid_rows = 'Setting',
                                                    steps = "confounder dimension",
                                                    methods = methods_depicted,
                                                    steps_y_base = -0.1, steps_y_height = 0.05,,
                                                    legend_breaks = methods_depicted,  # Specify the desired order
                                                    legend_labels = methods_depicted,  # Custom labels if needed
                                                    x_name = "Sample size", y_name = "Rejection rate",
                                                    spu_x_shift = 1,
                                                    colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                                    line_linetypes = c(1,2,3),
                                                    point_size = 4,
                                                    line_size = 1.5,
                                                    point_shapes = c(19,17,15),
                                                    steps_values_annotate = TRUE, steps_annotation_size = 8, steps_color='grey31',steps_annotation_color='grey31',
                                                    hline_intercept = c(0,0.05),
                                                    y_expand_add = c(0.1,0.15/2),
                                                    y_breaks = seq(0,1,0.2),
                                                    line_alpha =0.6,
                                                    point_alpha = 0.8,
                                                    legend_name = "DNCIT",
                                                    base_size = 28,
                                                    replace_labels = list(
                                                      Setting = c(#'No'='T1E',
                                                        'Yes'='Power')
                                                    ),
                                                    grid_labeller = labeller(#'No'='T1E',
                                                      'Yes'='Power'),
                                                    post_processing = list(
                                                      add_custom_theme = list(
                                                        axis.text.x = ggplot2::element_text(angle = -90,
                                                                                            vjust = 0.5)
                                                      )
                                                    ))
p_conf_dim_ci_mod <- p_conf_dim_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),  # Remove x-axis label from first plot
        axis.text.x = element_blank(),   # Remove x-axis ticks from first plot
        axis.ticks.x = element_blank())
p_conf_dim_no_ci_mod <- p_conf_dim_no_ci +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        legend.key.size = unit(2, "cm"),  # Increase overall size of the legend symbols
        legend.spacing.x = unit(0.5, "cm"), # Increase horizontal spacing between legend items
        legend.spacing.y = unit(0.5, "cm")) + # Adjust vertical spacing if needed
  guides(
    colour = guide_legend(nrow = 2, keywidth = 2, keyheight = 2),  # Increase width and height of the legend symbols
    shape = guide_legend(nrow = 2, keywidth = 2, keyheight = 2)    # Same for shape legend if needed
  )

combined_plot <- plot_grid(
  p_conf_dim_ci_mod,  # First plot
  p_conf_dim_no_ci_mod,  # Second plot with x-axis labels
  ncol = 1,  # Stack them vertically
  align = "v",  # Align the vertical axes
  rel_heights = c(1, 1.2)  # Adjust relative heights if needed
)
y_axis_label <- ggdraw() +
  draw_label("Rejection Rate", x = 0.5, y = 0.6, angle = 90, vjust = 0.5, size = 24) +
  theme(panel.background = element_rect(fill = "white", colour = NA), 
        plot.background = element_rect(fill = "white", colour = NA))
combined_plot <- plot_grid(
  y_axis_label,  # Y-axis label on the left
  combined_plot,  # The stacked plots
  ncol = 2,  # Arrange in two columns (label + plots)
  rel_widths = c(0.05, 1)  # Give space for the y-axis label
)
# Apply guides() to modify the legend appearance (e.g., number of rows, key size)
p_conf_dim_no_ci_legend <- p_conf_dim_legend +  # Customize the legend
  guides(colour = guide_legend(nrow = 3, keywidth = 2, keyheight=1.2, override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(size = 20))

# Use cowplot's get_legend() to extract the modified legend
legend <- get_legend(p_conf_dim_no_ci_legend)
p_conf_dim <- plot_grid(
  combined_plot,  # Stacked plots
  legend,  # Legend below the plots
  ncol = 1,  # Legend below, so keep 1 column
  rel_heights = c(1, 0.1)  # Adjust height ratios if needed
)
#ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_dim_main_text.png'), p_conf_dim, width = 18, height = 16, dpi = 300)
print(p_conf_dim)






### 3) Runtime
# only relevant files for runtime plot should be in folder!
folder_path <- "/sc/home/marco.simnacher/dncitPaper/Results/CI/runtime_cit/seeds_1_200"
all_files <- list.files(folder_path, full.names = TRUE)

## for sample sizes
folder_path_reject <- "/sc/home/marco.simnacher/dncitPaper/Results/CI/rejection_rates/seeds_1_200"
all_files_cit <- list.files(folder_path_reject, full.names = TRUE)
all_files_cit <- all_files_cit[setdiff(1:length(all_files_cit), grep('2_0_1_0|2_3_1_0|3_0_1_0|3_3_1_0|4_0_1_0|4_3_1_0|5_0_1_0|5_3_1_0', all_files_cit))]

cit_patterns <- "WALD|RCOT|kpc_graph|FCIT|CMIknn|comets_pcm"
cit_files <- grep(cit_patterns, all_files, value=TRUE)
files_conf_dim <- grep("squared", cit_files, value=TRUE)

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(145, 256, 350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
               confounder = c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15'))
design <- expand.grid(params)
design <- design %>%
  mutate("fastsurfer_fastsurfer-RCOT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-WALD" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-kpc_graph" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-FCIT" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-CMIknn" = rep(NA, nrow(design)),
         "fastsurfer_fastsurfer-comets_pcm" = rep(NA, nrow(design)),
         "freesurfer-RCOT" = rep(NA, nrow(design)),
         "freesurfer-WALD" = rep(NA, nrow(design)),
         "freesurfer-kpc_graph" = rep(NA, nrow(design)),
         "freesurfer-FCIT" = rep(NA, nrow(design)),
         "freesurfer-CMIknn" = rep(NA, nrow(design)),
         "freesurfer-comets_pcm" = rep(NA, nrow(design)),
         "condVAE-RCOT" = rep(NA, nrow(design)),
         "condVAE-WALD" = rep(NA, nrow(design)),
         "condVAE-kpc_graph" = rep(NA, nrow(design)),
         "condVAE-FCIT" = rep(NA, nrow(design)),
         "condVAE-CMIknn" = rep(NA, nrow(design)),
         "condVAE-comets_pcm" = rep(NA, nrow(design)),
         "medicalnet-RCOT" = rep(NA, nrow(design)),
         "medicalnet-WALD" = rep(NA, nrow(design)),
         "medicalnet-kpc_graph" = rep(NA, nrow(design)),
         "medicalnet-FCIT" = rep(NA, nrow(design)),
         "medicalnet-CMIknn" = rep(NA, nrow(design)),
         "medicalnet-comets_pcm" = rep(NA, nrow(design)))

embedding_maps <- c('fastsurfer_fastsurfer', 'freesurfer', 'condVAE', 'medicalnet')
dncits <- c('RCOT', 'WALD', 'kpc_graph', 'FCIT', 'CMIknn', 'comets_pcm')
confounder <- c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15')
for(dncit in dncits){
  for (embedding in embedding_maps){
    for (conf in confounder){
      if(conf == 'ukb_z1_'){
        files_dncit <- grep(dncit, all_files, value = TRUE)
        files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
        files <- grep(conf, files_dncit_conf, value=TRUE)
        files <- grep("squared", files, value=TRUE)
        df_runtimes <- read.csv(files, header = TRUE, sep = ",")

        df_runtimes <- colMeans(df_runtimes)

        files_dncit <- grep(dncit, all_files_cit, value = TRUE)
        files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
        #files <- grep(conf, files_dncit_conf, value=TRUE)
        df_cit <- read.csv(files_dncit_conf[1], header = TRUE, sep = ",")
        log_runtimes <- log(df_runtimes[2:length(df_runtimes)])
        df_cit[,2] <- log_runtimes


        col <- grepl(embedding, colnames(design)) & grepl(dncit, colnames(design))
        design[design$confounder==conf & design$sample_sizes %in% df_cit[,1], col] <- df_cit[,2]
      }else{
        if (TRUE){
          files_dncit <- grep(dncit, all_files, value = TRUE)
          files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
          files <- grep(conf, files_dncit_conf, value=TRUE)
          files <- grep("squared", files, value=TRUE)
          df_runtimes <- read.csv(files, header = TRUE, sep = ",")
          
          df_runtimes <- colMeans(df_runtimes)

          files_dncit <- grep(dncit, all_files_cit, value = TRUE)
          files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
          #files <- grep(conf, files_dncit_conf, value=TRUE)
          df_cit <- read.csv(files_dncit_conf[1], header = TRUE, sep = ",")
          log_runtimes <- log(df_runtimes[2:length(df_runtimes)])
          df_cit[,2] <- log_runtimes

          col <- grepl(embedding, colnames(design)) & grepl(dncit, colnames(design))
          design[design$confounder==conf & design$sample_sizes %in% df_cit[,1], col] <- df_cit[,2]
        }
      }

    }
  }
}
design_runtime <- design

##nested loop plot

design_runtime$confounder <- rep(c(1,2,4,6,10,15), each=10)
colnames(design_runtime) <- c("sample_sizes", "confounder dimension",
                      "Fastsurfer-RCOT", "Fastsurfer-WALD", "Fastsurfer-CPT_KPC", "Fastsurfer-FCIT", "Fastsurfer-CMIknn", "Fastsurfer-PCM",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "Freesurfer-CPT_KPC", "Freesurfer-FCIT", "Freesurfer-CMIknn", "Freesurfer-PCM",
                      "cVAE-RCOT", "cVAE-WALD", "cVAE-CPT_KPC", "cVAE-FCIT", "cVAE-CMIknn", "cVAE-PCM",
                      "MedicalNet-RCOT", "MedicalNet-WALD", "MedicalNet-CPT_KPC", "MedicalNet-FCIT", "MedicalNet-CMIknn", "MedicalNet-PCM")
custom_order <- c("sample_sizes", "confounder dimension",
                  "Fastsurfer-RCOT", "Freesurfer-RCOT","cVAE-RCOT", "MedicalNet-RCOT",
                  "Fastsurfer-CPT_KPC", "Freesurfer-CPT_KPC", "cVAE-CPT_KPC", "MedicalNet-CPT_KPC",
                  "Fastsurfer-FCIT", "Freesurfer-FCIT", "cVAE-FCIT", "MedicalNet-FCIT",
                  "Fastsurfer-CMIknn","Freesurfer-CMIknn", "cVAE-CMIknn", "MedicalNet-CMIknn",
                  "Fastsurfer-PCM", "Freesurfer-PCM","cVAE-PCM", "MedicalNet-PCM",
                  "Fastsurfer-WALD", "Freesurfer-WALD","cVAE-WALD", "MedicalNet-WALD")
#resort columns
design_runtime <- design_runtime[, custom_order]
design_runtime <- design_runtime[design_runtime$'confounder' %in% c(1,10),]
design_runtime <- design_runtime %>% select(-contains("Fastsurfer"))
# without PCM
#design_runtime <- design_runtime %>% select(-contains("PCM"))
methods_depicted <- colnames(design_runtime)[-c(1:2)]
y_lims = c(min(design_runtime[,-c(1,2)], na.rm = TRUE)-0.4, max(design_runtime[,-c(1,2)], na.rm = TRUE)+0.2)
p_conf_dim_runtime <- looplot::nested_loop_plot(resdf = design_runtime,
                                       x = "sample_sizes",
                                       steps  = "confounder dimension",
                                       methods = methods_depicted,
                                       legend_breaks = methods_depicted,  # Specify the desired order
                                       legend_labels = methods_depicted,  # Custom labels if needed
                                       steps_y_base = y_lims[1], steps_y_height = 0.15,
                                       x_name = "Sample size", y_name = expression("log"[10] * "(Runtime in s)"),
                                       spu_x_shift = 1,
                                       colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                       line_linetypes = rep(c(1,2,3), 5),
                                       point_size = 4,
                                       line_size = 1.5,
                                       point_shapes = rep(c(19,17,15),5),
                                       steps_values_annotate = TRUE, steps_annotation_size = 8, steps_color='grey31',steps_annotation_color='grey31',
                                       hline_intercept = y_lims[1]+0.4,
                                       y_expand_add = c(0.4,0.15),
                                       line_alpha =0.6,
                                       point_alpha = 0.8,
                                       ylim = y_lims,
                                       na_rm = FALSE,
                                       legend_name = "DNCIT",
                                       base_size = 38,
                                       post_processing = list(
                                         add_custom_theme = list(
                                           axis.text.x = ggplot2::element_text(angle = -90,
                                                                               vjust = 0.5)
                                         )
                                       ))
p_conf_dim_legend <- looplot::nested_loop_plot(resdf = design_runtime,
                                               x = "sample_sizes",
                                               grid_rows = 'Setting',
                                               steps = "confounder dimension",
                                               methods = methods_depicted,
                                               steps_y_base = -0.1, steps_y_height = 0.05,,
                                               legend_breaks = methods_depicted,  # Specify the desired order
                                               legend_labels = methods_depicted,  # Custom labels if needed
                                               x_name = "Sample size", y_name = "Rejection rate",
                                               spu_x_shift = 1,
                                               colors = palet_discrete[rep(c(1,2,7,4,3,6), each=3)],
                                               line_linetypes = c(1,2,3),
                                               point_size = 4,
                                               line_size = 1.5,
                                               point_shapes = c(19,17,15),
                                               steps_values_annotate = TRUE, steps_annotation_size = 8, steps_color='grey31',steps_annotation_color='grey31',
                                               hline_intercept = c(0,0.05),
                                               y_expand_add = c(0.1,0.15/2),
                                               y_breaks = seq(0,1,0.2),
                                               line_alpha =0.6,
                                               point_alpha = 0.8,
                                               legend_name = "DNCIT",
                                               base_size = 31,
                                               replace_labels = list(
                                                 Setting = c(#'No'='T1E',
                                                   'Yes'='Power')
                                               ),
                                               grid_labeller = labeller(#'No'='T1E',
                                                 'Yes'='Power'),
                                               post_processing = list(
                                                 add_custom_theme = list(
                                                   axis.text.x = ggplot2::element_text(angle = -90,
                                                                                       vjust = 0.5)
                                                 )
                                               ))

# Apply guides() to modify the legend appearance (e.g., number of rows, key size)
p_conf_dim_no_ci_legend <- p_conf_dim_legend +  # Customize the legend
  guides(colour = guide_legend(nrow = 3, keywidth = 2, keyheight=1.2, override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 19),
        legend.title = element_text(size = 20))
# Use cowplot's get_legend() to extract the modified legend
legend <- get_legend(p_conf_dim_no_ci_legend)
p_conf_dim_runtime <- p_conf_dim_runtime + theme(legend.position = 'none')
p_conf_dim_runtime <- plot_grid(
  p_conf_dim_runtime,  # Stacked plots
  legend,  # Legend below the plots
  ncol = 1,  # Legend below, so keep 1 column
  rel_heights = c(1, 0.1)  # Adjust height ratios if needed
)
ggplot2::ggsave(paste0(path_to_save_nested_loop_plots, '/all_dncits_nested_loop_conf_dim_runtime.png'), p_conf_dim_runtime, width = 19, height = 16, dpi = 300)
print(p_conf_dim_runtime)



#####QQ plots
palet_discrete <- paletteer::paletteer_d("ggthemes::Classic_10_Medium")
## conf dim 6, conf relationship complex
folder_path_qq <- "/sc/home/marco.simnacher/dncitPaper/Results/No_CI/p-values/seeds_1_200"
#folder_path_qq <- "Results\\No_CI\\p-values"
all_files_qq <- list.files(folder_path_qq, full.names = TRUE)
all_files_qq <- all_files_qq[setdiff(1:length(all_files_qq), grep('1_0_0.5_0', all_files_qq))]

### 1) Data preparation for nested loop over
# loop 1: confounder dimension (1,2,4,6,10,15) [APPENDIX] vs confounder dimension (1,2,10) [MAIN TEXT]
# loop 2: sample size
# fixed confounder relationship (squared terms of all continuous confounders)
cit_patterns <- "WALD|RCOT|kpc_graph|FCIT|CMIknn|comets_pcm"
cit_files <- grep(cit_patterns, all_files_qq, value=TRUE)
files_qq_squared <- grep("squared", cit_files, value=TRUE)
files_qq_dim_6 <- grep("ukb_z1_", files_qq_squared, value=TRUE)
files_qq <- grep("freesurfer", files_qq_dim_6, value=TRUE)
#files_qq_squared <- grep("squared", cit_files, value=TRUE)
#files_qq_dim_6 <- grep("ukb_z1_", files_qq_squared, value=TRUE)
#files_qq <- grep("freesurfer", files_qq_dim_6, value=TRUE)

dncits <- c('RCOT', 'kpc_graph', 'FCIT', 'CMIknn', 'WALD', 'comets_pcm')
p_values <- list()
for(dncit in dncits){
    files_dncit <- grep(dncit, files_qq, value = TRUE)
    df <- read.csv(files_dncit, header = TRUE, sep = ",")
    df <- df[,-c(4,6,8,9)]
    p_values[[dncit]] <- df
}

sample_sizes <- c(145, 256, 460, 1100, 5000, 10000)


# Determine the maximum number of columns (excluding the first column) across all data frames
max_cols <- max(sapply(p_values, function(df) ncol(df) - 1))
# Calculate the number of rows (number of dncits) and columns (maximum columns in data frames)
n_rows <- length(p_values)
n_cols <- max_cols
# List to store all the plots
plot_list <- list()

# Generate a palette with enough colors
colors <- setNames(palet_discrete[c(1,2,7,4,3,6)], dncits)

# Loop over each data frame in p_values
for (dncit_name in names(p_values)) {
  df <- p_values[[dncit_name]]
  num_cols <- ncol(df) - 1  # Exclude the first column
  plots_row <- list()       # Store plots for the current row

  # Loop over the maximum number of columns
  for (col_idx in 2:(max_cols + 1)) {
    if (dncit_name == "WALD"){
      if (col_idx == 2){
        # Create an empty plot
        p <- ggplot() +
          ggtitle(paste(dncit_name, "- Column", col_idx - 1)) +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "Theoretically inapplicable", size = 5, hjust = 0.5)
        p <- p + theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_blank(),
          text = element_text(size = 14)
        )
        p <- p + theme(
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          text = element_text(size = 14)
        )

        plots_row[[length(plots_row) + 1]] <- p
        # Create an empty plot
        p <- ggplot() +
          ggtitle(paste(dncit_name, "- Column", col_idx)) +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "Theoretically inapplicable", size = 5, hjust = 0.5)
        p <- p + theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_blank(),
          text = element_text(size = 14)
        )
        p <- p + theme(
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          text = element_text(size = 14)
        )

        plots_row[[length(plots_row) + 1]] <- p

        # Extract p-values from the current column
        p_vals <- df[[col_idx]]
        # Remove NA values
        p_vals <- p_vals[!is.na(p_vals)]
        ks <- round(ks.test(p_vals, "punif")[1]$statistic,3)

        # Create a data frame for plotting
        plot_df <- data.frame(
          Observed = sort(p_vals),
          Theoretical = qunif(ppoints(length(p_vals)))
        )

        # Create the QQ plot
        p <- ggplot(plot_df, aes(sample = Observed)) +
          stat_qq(distribution = stats::qunif, color = colors[dncit_name]) +
          geom_abline(slope = 1, intercept = 0, color = "black") +
          ggtitle(paste(dncit_name, "- Column", col_idx + 1)) +
          xlab(ks) +
          ylab("Sample Quantiles") +
          theme_minimal() +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
        p <- p + theme(
          #axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_blank(),
          text = element_text(size = 14)
        )

        plots_row[[length(plots_row) + 1]] <- p
      }else{
        if (col_idx <= ncol(df)) {
          # Extract p-values from the current column
          p_vals <- df[[col_idx]]
          # Remove NA values
          p_vals <- p_vals[!is.na(p_vals)]
          ks <- round(ks.test(p_vals, "punif")[1]$statistic,3)

          # Create a data frame for plotting
          plot_df <- data.frame(
            Observed = sort(p_vals),
            Theoretical = qunif(ppoints(length(p_vals)))
          )

          # Create the QQ plot
          p <- ggplot(plot_df, aes(sample = Observed)) +
            stat_qq(distribution = stats::qunif, color = colors[dncit_name]) +
            geom_abline(slope = 1, intercept = 0, color = "black") +
            ggtitle(paste(dncit_name, "- Column", col_idx + 1)) +
            xlab(ks) +
            ylab("Sample Quantiles") +
            theme_minimal() +
            coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
          p <- p + theme(
            axis.title.y = element_blank(),
            #axis.title.x = element_blank(),
            plot.title = element_blank(),
            text = element_text(size = 14)
          )
          if(col_idx!=ncol(df)){
            p <- p + theme(
              #axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_blank(),
              text = element_text(size = 14)
            )
          }else{
            p <- p + theme(
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_blank(),
              text = element_text(size = 14)
            )
          }
          plots_row[[length(plots_row) + 1]] <- p
        }
      }
    }else{
      if (col_idx <= ncol(df)) {
        # Extract p-values from the current column
        p_vals <- df[[col_idx]]
        # Remove NA values
        p_vals <- p_vals[!is.na(p_vals)]
        ks <- round(ks.test(p_vals, "punif")[1]$statistic,3)

        # Create a data frame for plotting
        plot_df <- data.frame(
          Observed = sort(p_vals),
          Theoretical = qunif(ppoints(length(p_vals)))
        )

        # Create the QQ plot
        p <- ggplot(plot_df, aes(sample = Observed)) +
          stat_qq(distribution = stats::qunif, color = colors[dncit_name]) +
          geom_abline(slope = 1, intercept = 0, color = "black") +
          ggtitle(paste(dncit_name, "- Column", col_idx - 1)) +
          xlab(ks) +
          ylab(sample_sizes[col_idx-1]) +
          theme_minimal()+
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
        p <- p + theme(
          axis.title.y = element_blank(),
          #axis.title.x = element_blank(),
          plot.title = element_blank(),
          text = element_text(size = 14),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA)
        )
        if(col_idx!=(max_cols+1) & dncit_name!='RCOT'){
          p <- p + theme(
            #axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_blank(),
            text = element_text(size = 14)
          )
        }else if(col_idx == (max_cols+1) & dncit_name!='RCOT'){
            p <- p + theme(
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_blank(),
              text = element_text(size = 14)
            )
        }else if(col_idx != (max_cols+1) & dncit_name=='RCOT'){
          p <- p + theme(
            #axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_blank(),
            text = element_text(size = 14)
          )
        }
        if(dncit_name=='RCOT'){
          p <- p + theme(
            axis.title.y = element_text(angle=90,margin = margin(t = 0, r = 5, b = 0, l = 0))
          )
        }
        plots_row[[length(plots_row) + 1]] <- p
      } else {
        # Create an empty plot
        p <- ggplot() +
          ggtitle(paste(dncit_name, "- Column", col_idx - 1)) +
          theme_void() +
          annotate("text", x = 0.5, y = 0.5, label = "Computational restrictions", size = 5, hjust = 0.5)
        p <- p + theme(
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          text = element_text(size = 14),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA)
        )
        if(col_idx!=max_cols){
          p <- p + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_blank(),
            text = element_text(size = 14)
          )
        }else{
          p <- p + theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_blank(),
            text = element_text(size = 14)
          )
        }

        plots_row[[length(plots_row) + 1]] <- p
      }
  }
  # Add the row of plots to the main plot list
  }
  plot_list <- c(plot_list, plots_row)
}

col_names <- c('Freesurfer-RCOT', 'Freesurfer-CPT-KPC', 'Freesurfer-FCIT', 'Freesurfer-CMIknn', 'Freesurfer-WALD', 'Freesurfer-PCM')
# Create text grobs for the titles
title_grobs <- lapply(col_names, function(title) {
  textGrob(title, gp = gpar(fontsize = 18), just = "center")
})

# Combine the title grobs and plot list
grob_list <- c(title_grobs, plot_list)

# Create a layout matrix to arrange plots row-wise
layout_matrix <- rbind(matrix(1:n_rows, nrow=1, ncol=n_rows),matrix(c((n_rows + 1):(n_rows + n_rows * n_cols)), nrow=n_cols, byrow=FALSE))


# Arrange all the plots into a grid using the layout matrix
grid_plots <- gridExtra::grid.arrange(
  grobs = grob_list,
  layout_matrix = layout_matrix,
  bottom = "Theoretical Quantiles",
  left = "Sample Quantiles",
  heights = unit.c(unit(1, "lines"), unit(1.1,"null"), rep(unit(1, "null"), n_cols-1))
)
# Display the grid of plots
grid::grid.newpage()
grid::grid.draw(grid_plots)

# Save as PNG
#ggsave(paste0(path_to_save_nested_loop_plots,"/qq_plot_dim6_squared_freesurfer_ci.png"), grid_plots, width = 12, height = 12, dpi = 300)
# Save as PDF
ggsave(paste0(path_to_save_nested_loop_plots,"/qq_plot_dim1_squared_freesurfer_ci.pdf"), grid_plots, width = 14, height = 14)
# ggsave(paste0(path_to_save_nested_loop_plots,"qq_plot_dim1_squared_freesurfer_no_ci.png"), grid_plots, width = 12, height = 12)
ggsave(paste0(path_to_save_nested_loop_plots,"/qq_plot_dim1_squared_freesurfer_no_ci.pdf"), grid_plots, width = 14, height = 14)

##### Tables for detailed results
library(xtable)
rownames(design) <- NULL
xtable(design[,-c(ncol(design))], align=c("lll|rrrrrrrrrrrrrrrrrr"))
heatmap(as.matrix(design[,-c(1,2,ncol(design))]))
