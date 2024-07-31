### Definition of simulation settings and run over all settings for all DNCITs
setwd('inst')
# Settings
#### args =['/CI/ or /No_CI/', post_non_lin, eps_sigmaX, eps_sigmaY, eps_sigmaZ, embedding_orig, embedding_obs, confounder, g_z, CIT, CIT_params]
#c("/CI/", "1", "0", "1", "0", "fastsurfer", "fastsurfer", "AS", "linear", "RCOT", "1")
settings <- data.frame(
  dependence = rep(c('/CI/', '/No_CI/'), each = 120),
  fct_relation = rep(rep(c('1','2','3','4','5'), each=4), 12),
  eps_sigmaX = rep(rep(c(0, 3, 0, 0), 10), 6),
  eps_sigmaY = rep(1, 80),
  eps_sigmaZ = rep(0, 80),
  embedding_orig = rep('fastsurfer', 80),
  embedding_obs = rep(c('fastsurfer', 'noisy', 'freesurfer', 'condVAE'), 60),
  confounder = rep(rep(c("ukb_z1","ukb_z2", "ukb_z4", "ukb_z6", "ukb_z10", "ukb_z15"), each=20),2),
  response = rep(c('linear', 'squared', 'realistic', 'breakpoint3'), 60)
)

cits <- c('WALD', 'RCOT 1')#, 'kpc_graph 2 10')#'CMIknn')#, 'FCIT'
settings <- settings[1:4,]
# Run 'sim_ukb_brainmri.R' for all settings
for (cit in cits){
  for (i in 1:nrow(settings)) {
    setting <- paste(settings[i,], collapse = " ")
    command <- paste("Rscript raw_code/sim_ukb_brainmri.R", setting, cit)

    # Execute the command
    system(command)
  }
}

