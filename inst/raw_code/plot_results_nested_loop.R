##### Plot rejection rates of CITs in nested loop plots (cf. Rücker et al. 2014)
##color palettes for plotting
palet_discrete <- paletteer::paletteer_d("colorBlindness::Blue2Orange10Steps")


#### CI and No CI
folder_path_ci <- "M:\\CITs\\Application\\UKB_data\\Results\\CI\\rejection_rates"
folder_path_no_ci <- "M:\\CITs\\Application\\UKB_data\\Results\\No_CI\\rejection_rates"
all_files_ci <- list.files(folder_path_ci, full.names = TRUE)
all_files_no_ci <- list.files(folder_path_no_ci, full.names = TRUE)
all_files <- union(all_files_ci, all_files_no_ci)

### 1) nested loop over
# loop 1: confounder dimension (1,2,4,6,10,15)
# loop 2: sample size
# fixed confounder relationship (squared terms of all continuous confounders)
# plot Deep-WALD and Deep-RCOT for all four embedding maps (freesurfer, fastsurfer, condVAE, noisy fastsurfer)
squared_conf_files <- all_files[grep("squared", all_files)]
wald_files <- all_files[grep("WALD", all_files)]
rcot_files <- all_files[grep("RCOT", all_files)]
files_conf_dim <- union(intersect(squared_conf_files, wald_files), intersect(squared_conf_files, rcot_files))

# Result tab for each DNCIT
params <- list(sample_sizes = as.factor(c(350, 460, 825, 1100, 1475, 1964, 5000, 10000)),
      confounder = c('ukb_z1_', 'ukb_z2', 'ukb_z4', 'ukb_z6', 'ukb_z10', 'ukb_z15'),
      dependence = c('CI','No_CI'))
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
for(dependence in c('Results\\\\CI','No_CI')){
  for(dncit in dncits){
    for (embedding in embedding_maps){
      for (conf in confounder){
        files_ci <- grep(dependence, files_conf_dim, value = TRUE)
        files_dncit <- grep(dncit, files_ci, value = TRUE)
        files_dncit_conf <- grep(embedding, files_dncit, value = TRUE)
        files <- grep(conf, files_dncit_conf, value=TRUE)
        if(length(files)>1){
          files <- grep('1_0_1_0', files, value=TRUE)
        }
        df <- read.csv(files, header = TRUE, sep = ",")
        col <- grepl(embedding, colnames(design)) & grepl(dncit, colnames(design))
        design[design$confounder==conf, col] <- df[,2]
      }
    }
  }
}


design$confounder <- rep(c(1,2,4,6,10,15), each=8)
colnames(design) <- c("sample_sizes", "confounder dimension", "Fastsurfer-RCOT", "Fastsurfer-WALD",
                      "Freesurfer-RCOT", "Freesurfer-WALD", "condVAE-RCOT", "condVAE-WALD")

## nested loop plot
p = looplot::nested_loop_plot(resdf = design,
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
print(p)
