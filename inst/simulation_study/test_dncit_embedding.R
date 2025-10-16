##### Testing for simulation  
library(devtools)
load_all()

# Load embeddings using reticulate (numpy only - should work fine)
library(reticulate)
np <- import("numpy")

######### Data checking #########
# load IXI image_path.csv after running prepare_ixi_from_tar_xls.py
readRenviron(".Renviron")
Sys.getenv(c("UKB_PATH","IXI_PATH"))
ukb_path <- Sys.getenv("UKB_PATH", unset = NA)
script_dir <- "/home/RDC/simnacma/Coding/dncitPaper/inst/learn_embedding"

ukb_data <- read.csv(paste0(ukb_path, "t1_paths.csv"))
### Example brain MRI plot
library(RNifti)  
idx <- 1
nii_path <- ixi_data$path[idx]
nii <- RNifti::readNifti(nii_path)
# Choose a mid axial slice
d   <- dim(nii)
z   <- ceiling(d[3]/2)
# Output file name
id <- ixi_data$id[idx]
outfile <- file.path(ixi_path, sprintf("t1/example/%s_axial_z%03d.png", id, z))
# Extract slice, normalize for display, flip for conventional orientation, and save
sl <- nii[,,z, drop = TRUE]
sl <- sl - min(sl, na.rm = TRUE)
sl <- if (max(sl, na.rm = TRUE) > 0) sl / max(sl, na.rm = TRUE) else sl
png(outfile, width = 1200, height = 1200, res = 150)
par(mar = c(0,0,0,0), bg = "white")
image(t(apply(sl, 2, rev)), axes = FALSE, col = gray.colors(256), useRaster = TRUE, asp = 1)
dev.off()



# set up as in sim_run_settings.R and sim_ukb_brainmri.R
args <- c("/CI/", "1", "0", "1", "0", "fastsurfer", "fastsurfer", "ukb_z4", "linear", "ccit")
n_cits <- 1
cit <- c(args[10])
print(args)

n_seeds = 1:200

post_non_lin=as.numeric(args[2])
eps_sigmaX=as.numeric(args[3])
eps_sigmaY=as.numeric(args[4])
eps_sigmaZ=as.numeric(args[5])
embedding_orig=args[6]
embedding_obs=args[7]
confounder=args[8]                                        
g_z=args[9]
idx_beta2=NULL
idx_sample = 1
seed=1

## test data generation
set.seed(seed)

##### Test DNCIT-specific embedding maps
## prepare data
path_to_ukb_data <- Sys.getenv("UKB_PATH", unset = NA)

if (is.na(path_to_ukb_data) || path_to_ukb_data == "") {
  stop("❌ Environment variable UKB_PATH is not set. Please define it in your .Renviron.")
}
X_orig <- load_X_orig(path_to_ukb_data,embedding_orig)
X_obs <- load_X_obs(path_to_ukb_data,embedding_obs, embedding_orig, X_orig, eps_sigmaX)
Z <- load_Z(path_to_ukb_data,confounder)

# take only common IDs of X, Z
colnames(X_orig)[1] <- 'id'
colnames(X_obs)[1] <- 'id'
colnames(Z)[1] <- 'id'
X_orig <- as.data.frame(X_orig)
X_obs <- as.data.frame(X_obs)
Z <- as.data.frame(Z)

common_ids <- Reduce(intersect, list(X_orig$id, X_obs$id, Z$id))

subset_X_obs <- X_obs[X_obs$id %in% common_ids, ]
subset_X_orig <-  X_orig[X_orig$id %in% common_ids, ]
subset_Z <- Z[Z$id %in% common_ids, ]
#check if ids are equal
stopifnot(all.equal(subset_X_orig[,1], subset_Z[,1]))

#sample rows
idx <- sample(1:nrow(subset_Z), n_sample[[idx_sample]])
X_orig <- subset_X_orig[idx,-c(1)]
Z <- subset_Z[idx,-c(1), drop=FALSE]

epsZ <- matrix(stats::rnorm((nrow(Z)*ncol(Z)), 0, eps_sigmaZ), nrow=nrow(Z), ncol=ncol(Z))
Z <- Z+epsZ

  #remove zero columns and multicollinearity in one-hot-encoding of sites (due to subsampling)
  site_columns <- grep("^site", colnames(Z), value = TRUE)
  for (site in site_columns){
    site_sum <- sum(Z[[site]])
    if(site_sum == 0){
      Z <- Z[, !names(Z) %in% site]
    }
  }
  site_columns <- grep("^site", colnames(Z), value = TRUE)
  # Sum the 'site' columns row-wise
  site_sum <- sum(Z[site_columns])
  #remove one site column
  if(site_sum == nrow(Z)){
    for(site in site_columns){
      if(is_binary(Z[[site]])){
        Z <- Z[, !names(Z) %in% site]
        break
      }
    }
  }

  #Standardize
  # scale only continuous confounders
  Z <- Z %>% dplyr::mutate(dplyr::across(dplyr::where(function(x) !is_binary(x)), scale))

  X_orig <- scale(X_orig)
  X_obs <- scale(X_obs)

  if(is.null(beta2s)){
    Y <- y_from_xz(Z, eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z)
  }else{
    Y <- y_from_xz(Z, eps_sigmaY, X=X_orig, beta2s=beta2s, idx_beta2=idx_beta2, post_non_lin=post_non_lin, g_z=g_z)
  }

