### Definition of simulation settings and run over all settings for all DNCITs
setwd('inst')
# Settings
#### args =['/CI/ or /No_CI/', post_non_lin, eps_sigmaX, eps_sigmaY, eps_sigmaZ, embedding_orig, embedding_obs, confounder, g_z, CIT, CIT_params]
if(FALSE){
  args <- c("/CI/", "1", "0", "1", "0", "fastsurfer", "fastsurfer", "ukb_z4", "linear", "comets_gcm")
  idx_sample=idx_beta2=i=1
  n_sample = list(350, 460, 825, 1100, 1475, 1964, 5000, 10000)
  XYZ_list <- dncitPaper::data_gen(seed=i, idx_sample=idx_sample, n_sample=n_sample, idx_beta2=NULL, beta2s=NULL,
                                                                            post_non_lin=as.numeric(args[2]), eps_sigmaX=as.numeric(args[3]), eps_sigmaY=as.numeric(args[4]),
                                                                            eps_sigmaZ=as.numeric(args[5]), embedding_orig=args[6], embedding_obs=args[7],
                                                                            confounder=args[8], g_z=args[9])
  X <- as.matrix(XYZ_list[[1]])
  Y <- as.matrix(XYZ_list[[2]])
  Z <- as.matrix(XYZ_list[[3]])

  #### old, all settings
  settings <- data.frame(
    dependence = rep(c('/CI/', '/No_CI/'), each = 480),
    fct_relation = rep(rep(c('1','2','3','4','5'), each=4), 48),
    eps_sigmaX = rep(rep(c(0, 3, 0, 0), 10), 24),
    eps_sigmaY = rep(1, 960),
    eps_sigmaZ = rep(0, 960),
    embedding_orig = rep('fastsurfer', 960),
    embedding_obs = rep(c('fastsurfer', 'noisy', 'freesurfer', 'condVAE'), 240),
    confounder = rep(rep(c("ukb_z1","ukb_z2", "ukb_z4", "ukb_z6", "ukb_z10", "ukb_z15"), each=80),2),
    response = rep(rep(c('linear', 'squared', 'realistic', 'breakpoint3'), each=20), 12)
  )
}

cits <- c('RCOT 1')#'FCIT', 'WALD', 'cpi', CMIknn', 'comets_gcm', 'comets_pcm', 'pred_cit','kpc_graph 2 10'
#settings <- settings#[c(41, 61,81,101,121,141,161),]
#runtime settings
settings_dim_z <- data.frame(dependence = rep(c('/CI/', '/No_CI/'), each = 24),
                             fct_relation = rep('1', 48),
                             eps_sigmaX = rep(c(0, 3, 0, 0), 12),
                             eps_sigmaY = rep(1, 48),
                             eps_sigmaZ = rep(0, 48),
                             embedding_orig = rep('fastsurfer', 48),
                             embedding_obs = rep(c('fastsurfer', 'noisy', 'freesurfer', 'condVAE'), 12),
                             confounder = rep(rep(c("ukb_z1","ukb_z2", "ukb_z4", "ukb_z6", "ukb_z10", "ukb_z15"), each=4),2),
                             response =  rep('squared', 48)
)
settings_g_z <- data.frame(dependence = rep(c('/CI/', '/No_CI/'), each = 12),
                           fct_relation = rep('1', 24),
                           eps_sigmaX = rep(c(0, 3, 0, 0), 6),
                           eps_sigmaY = rep(1, 24),
                           eps_sigmaZ = rep(0, 24),
                           embedding_orig = rep('fastsurfer', 24),
                           embedding_obs = rep(c('fastsurfer', 'noisy', 'freesurfer', 'condVAE'), 6),
                           confounder = rep("ukb_z6", 24),
                           response =  rep(rep(c('linear', 'squared', 'realistic'), each=4), 2)
)
settings <- rbind(settings_g_z, settings_dim_z)
## remove overlapping settings for g_z and dim_z
settings <- settings[-c(37:40, 61:64),]
# settings for the runtime
settings <- settings[c(25:28, 31, 35, 7, 39, 43),]
# settings <- settings[-c(1:16),]
# settings <- settings[settings$embedding_obs=="noisy" & settings$dependence=="/CI/",]
# Run 'sim_ukb_brainmri.R' for all settings
for (cit in cits){
  for (i in 1:nrow(settings)) {
    setting <- paste(settings[i,], collapse = " ")
    command <- paste("Rscript raw_code/sim_ukb_brainmri.R", setting, cit)

    # Execute the command
    system(command)
  }
}


#### To obtain .txt files which can be copied into the .sh files
settings_g_z <- data.frame(cit = rep('Rscript'), script=rep('raw_code/sim_ukb_brainmri.R'),
  dependence = rep(c('/CI/', '/No_CI/'), each = 12),
  fct_relation = rep('1', 24),
  eps_sigmaX = rep(c(0, 3, 0, 0), 6),
  eps_sigmaY = rep(1, 24),
  eps_sigmaZ = rep(0, 24),
  embedding_orig = rep('fastsurfer', 24),
  embedding_obs = rep(c('fastsurfer', 'noisy', 'freesurfer', 'condVAE'), 6),
  confounder = rep("ukb_z6", 24),
  response =  rep(rep(c('linear', 'squared', 'realistic'), each=4), 2)
)

settings_dim_z <- data.frame(cit = rep('Rscript'), script=rep('raw_code/sim_ukb_brainmri.R'),
  dependence = rep(c('/CI/', '/No_CI/'), each = 24),
  fct_relation = rep('1', 48),
  eps_sigmaX = rep(c(0, 3, 0, 0), 12),
  eps_sigmaY = rep(1, 48),
  eps_sigmaZ = rep(0, 48),
  embedding_orig = rep('fastsurfer', 48),
  embedding_obs = rep(c('fastsurfer', 'noisy', 'freesurfer', 'condVAE'), 12),
  confounder = rep(rep(c("ukb_z1","ukb_z2", "ukb_z4", "ukb_z6", "ukb_z10", "ukb_z15"), each=4),2),
  response =  rep('squared', 48)
)

settings_final <- rbind(settings_g_z, settings_dim_z)
settings_final$CIT <- 'RCOT 1'
write.table(settings_final, file = "raw_code/slurm/settings_rcot.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

settings_final <- rbind(settings_g_z, settings_dim_z)
settings_final$CIT <- 'CMIknn'
write.table(settings_final, file = "raw_code/slurm/settings_cmiknn.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

settings_final <- rbind(settings_g_z, settings_dim_z)
settings_final$CIT <- 'kpc_graph 2 10'
write.table(settings_final, file = "raw_code/slurm/settings_kpc_graph 2 10.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

settings_final <- rbind(settings_g_z, settings_dim_z)
settings_final$CIT <- 'WALD'
write.table(settings_final, file = "raw_code/slurm/settings_wald.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)

settings_final <- rbind(settings_g_z, settings_dim_z)
settings_final$CIT <- 'FCIT'
write.table(settings_final, file = "raw_code/slurm/settings_fcit.txt", sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
