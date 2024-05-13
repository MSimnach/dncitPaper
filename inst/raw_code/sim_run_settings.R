### Definition of simulation settings and run over all settings for all DNCITs
setwd('inst')
# Settings
#### args =['/CI/ or /No_CI/', post_non_lin, eps_sigmaX, eps_sigmaY, eps_sigmaZ, embedding_orig, embedding_obs, confounder, response, CIT, CIT_params]
#c("/CI/", "1", "0", "1", "0", "fastsurfer", "fastsurfer", "AS", "simulated", "RCOT", "1")
settings <- data.frame(
  dependence = rep(c('/CI/', '/No_CI/'), each = 24),
  fct_relation = rep(rep(c('1','4','5'), each=4), 4),
  eps_sigmaX = rep(rep(c(0, 3,0, 0), 6), 2),
  eps_sigmaY = rep(1, 48),
  eps_sigmaZ = rep(0, 48),
  embedding_orig = rep('freesurfer', 48),
  embedding_obs = rep(c('fastsurfer', 'noisy', 'freesurfer', 'condVAE'), 12),
  confounder = rep(rep(c("AS", "genes10"), each=12),2),
  response = rep('simulated', 48)
)

cits <- c('RCOT 1', 'WALD')#, 'kpc_graph 2 10', 'FCIT', 'CMIknn')
settings <- settings
# Run 'sim_ukb_brainmri.R' for all settings
for (cit in cits){
  for (i in 1:nrow(settings)) {
    setting <- paste(settings[i,], collapse = " ")
    command <- paste("Rscript raw_code/sim_ukb_brainmri.R", setting, cit)

    # Execute the command
    system(command)
  }
}

