run_python_safe <- function(script_path, args = character(), python = NULL) {
  stopifnot(file.exists(script_path))
  if (length(args) == 1) args <- strsplit(args, " +")[[1]]  # optional: allow single string

  # Prefer current conda env's python
  if (is.null(python) || !nzchar(python)) {
    cp <- Sys.getenv("CONDA_PREFIX", unset = "")
    if (nzchar(cp) && file.exists(file.path(cp, "bin", "python"))) {
      python <- file.path(cp, "bin", "python")
    } else {
      python <- Sys.which("python"); if (!nzchar(python)) python <- "python"
    }
  }

  # In case a previous helper script polluted the session:
  if (nzchar(Sys.getenv("LD_PRELOAD", ""))) Sys.unsetenv("LD_PRELOAD")

  out <- system2(python, c(normalizePath(script_path), args), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status"); if (is.null(status)) status <- 0L
  if (status != 0L) stop(sprintf("Python exited with status %s\n%s", status, paste(out, collapse = "\n")))
  out
}
##### Testing for simulation  
library(devtools)
load_all()
library(arrow)

# Load embeddings using reticulate (numpy only - should work fine)
#library(reticulate)
#np <- import("numpy")

######### Data checking #########
# load UKB image_path.csv 
readRenviron(".Renviron")
Sys.getenv(c("UKB_PATH","IXI_PATH"))
ukb_path <- Sys.getenv("UKB_PATH", unset = NA)
script_dir <- "inst/learn_embedding"

ukb_data <- read.csv(paste0(ukb_path, "/t1_paths_cit.csv"))
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
args <- c("/No_CI/", "1", "0", "0.05", "0", "fastsurfer", "scratch", "ukb_z4", "squared", "RCOT", "1")
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
beta2s=list(1)
idx_beta2=1
idx_sample = 1
seed=1
n_sample = list(550)

## test data generation
set.seed(seed)

##### Test DNCIT-specific embedding maps
## prepare data
path_to_ukb_data <- Sys.getenv("UKB_PATH", unset = NA)

# data loading functions
load_X_orig <- function(path_to_ukb_data, embedding_orig){
  if(embedding_orig == 'fastsurfer'){
    X <- data.table::fread(paste0(path_to_ukb_data, '/ukb_fastsurfer.csv'), header=TRUE, nThread = 1)
  }else if(embedding_orig == 'condVAE'){
    X <- data.table::fread(paste0(path_to_ukb_data, "/ukb_condVAE.csv"), header=TRUE, nThread = 1)
  }else if(embedding_orig == 'latentDiffusion'){
    X <- data.table::fread(paste0(path_to_ukb_data, "/ukb_latentDiffusion.csv"), header=TRUE, nThread = 1)
  }else if(embedding_orig == 'freesurfer'){
    X <- data.table::fread(paste0(path_to_ukb_data, "/ukb_freesurfer.csv"), header=TRUE, nThread = 1)
  }
  return(X)
}

load_X_obs <- function(path_to_ukb_data,embedding_obs, embedding_orig, X_orig, eps_sigmaX){
  if(embedding_obs==embedding_orig){
    X_obs <- X_orig
  }else if(embedding_obs == 'fastsurfer'){
    X_obs <- data.table::fread(paste0(path_to_ukb_data, "/ukb_fastsurfer.csv"), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'condVAE'){
    X_obs <- data.table::fread(paste0(path_to_ukb_data, "/ukb_condVAE.csv"), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'latentDiffusion'){
    X_obs <- data.table::fread(paste0(path_to_ukb_data, "/ukb_latentDiffusion.csv"), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'freesurfer'){
    X_obs <- data.table::fread(paste0(path_to_ukb_data, "/ukb_freesurfer.csv"), header=TRUE, nThread = 1)
  }else if(grepl('noisy',embedding_obs, fixed=TRUE)){
    epsX <- stats::rnorm(nrow(X_orig)*(ncol(X_orig)-1), 0,eps_sigmaX)
    X_obs <- cbind(X_orig[,1], scale(X_orig[,2:ncol(X_orig)])+epsX)
  }else if(embedding_obs %in% c('scratch', 'medicalnet_ft', 'medicalnet_ft_frozen')){
    X_obs <- data.table::fread(paste0(path_to_ukb_data, '/t1_paths_cit.csv'), header=TRUE, nThread = 1)
  }else if (embedding_obs == 'medicalnet'){
    X_obs_emb <- arrow::read_parquet(paste0(path_to_ukb_data, '/medicalnet_embeddings/embeddings.parquet'))
    X_obs <- as.data.frame(X_obs_emb)
  }
  return(X_obs)
}

load_Z <- function(path_to_ukb_data,confounder){
  if(confounder=='AS'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_Z.csv"), header=TRUE, nThread = 1)
  }else if(confounder == 'genes10'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_Z_genes10.csv"), header=TRUE, nThread = 1)
  }else if(confounder == 'ukb_z1'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z1_age.csv"), header=TRUE, nThread = 1)
  }else if(confounder == 'ukb_z2'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z2_agesex.csv"), header=TRUE, nThread = 1)
  }else if(confounder == 'ukb_z4'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z4_agesexsitesize.csv"), header=TRUE, nThread = 1)
  }else if(confounder == 'ukb_z6'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z6_agesexsitesizedateqc.csv"), header=TRUE, nThread = 1)
  }else if(confounder == 'ukb_z10'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z10_agesexsitesizedateqclocation.csv"), header=TRUE, nThread = 1)
  }else if(confounder == 'ukb_z15'){
    Z <- data.table::fread(paste0(path_to_ukb_data, "/ukb_z15_agesexsitesizedateqcgenes.csv"), header=TRUE, nThread = 1)
  }
}

# Function to check if a column has only 2 values (binary after sd e.g.)
is_binary <- function(x) {
  unique_values <- unique(x)
  length(unique_values) == 2
}

Z <- load_Z(path_to_ukb_data,confounder)
X_orig <- load_X_orig(path_to_ukb_data,embedding_orig)

# take only common IDs of X, Z
colnames(X_orig)[1] <- 'id'
colnames(Z)[1] <- 'id'
X_orig <- as.data.frame(X_orig)
Z <- as.data.frame(Z)

#Z <- Z[!duplicated(Z$id), ]
cit_ids <- data.table::fread(paste0(path_to_ukb_data, '/t1_paths_cit.csv'), header=TRUE, nThread = 1)$id
common_ids <- Reduce(intersect, list(X_orig$id, Z$id, cit_ids))

subset_X_orig <-  X_orig[match(common_ids, X_orig$id), ]
subset_Z <- Z[match(common_ids, Z$id), ]
#check if ids are equal
stopifnot(all.equal(subset_X_orig[,1], subset_Z[,1]))

#sample rows
idx <- sample(1:nrow(subset_Z), n_sample[[idx_sample]])
X_orig <- subset_X_orig[idx,]
Z <- subset_Z[idx,]

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
Z[,c(-1)] <- Z[,c(-1)] %>% dplyr::mutate(dplyr::across(dplyr::where(function(x) !is_binary(x)), scale))
# scale X_orig
X_orig[,c(-1)] <- scale(X_orig[,c(-1)])

if(is.null(beta2s)){
  Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z)
}else{
  Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, X=as.matrix(X_orig[,c(-1)]), beta2s=beta2s, idx_beta2=idx_beta2, post_non_lin=post_non_lin, g_z=g_z)
}

y_dir <- paste0(path_to_ukb_data, '/', args[1], '/', n_sample[[idx_sample]], '/', seed)
dir.create(y_dir, recursive = TRUE)
Y_id <- data.frame(id = X_orig$id, Y = Y[,1])
write.csv(Y_id, file.path(y_dir, 'Y.csv'), row.names = FALSE)


if(embedding_obs %in% c('fastsurfer', 'condVAE', 'latentDiffusion', 'freesurfer')){
  X_obs[,c(-1)] <- scale(X_obs[,c(-1)])
} else{
  # Reorder X_obs to match X_orig,Z,Y and assign Y to X_obs
  X_obs <- data.table::fread(paste0(path_to_ukb_data, '/t1_paths_cit.csv'), header=TRUE, nThread = 1)
  X_obs <- X_obs[match(X_orig$id, X_obs$id), ]
  stopifnot(all.equal(X_obs$id, Y_id$id))
  X_obs$y <- Y[,1]

  # Setup paths for train results and embeddings
  embedding_dir <- paste0(y_dir, '/', embedding_obs)
  embeddings_path <- file.path(embedding_dir, '/test_embeddings.parquet')
  index_path <- file.path(embedding_dir, '/test_embeddings_index.csv')
  preds_path <- file.path(embedding_dir, '/test_predictions.parquet')

  ## Train or skip if embeddings already exist
  if (dir.exists(embedding_dir)){
    cat("✅ Embedding directory already exists. Skipping training...\n")
  } else {
    dir.create(embedding_dir, recursive = TRUE)
    # save Y and data gen config
    # Save configuration
    config <- list(
      args = args,
      post_non_lin = post_non_lin,
      eps_sigmaX = eps_sigmaX,
      eps_sigmaY = eps_sigmaY,
      eps_sigmaZ = eps_sigmaZ,
      embedding_orig = embedding_orig,
      embedding_obs = embedding_obs,
      confounder = confounder,
      g_z = g_z,
      cit = cit,
      seed = seed,
      n_sample = n_sample[[idx_sample]],
      idx_sample = idx_sample,
      timestamp = Sys.time()
    )
    # Save configuration as JSON for easy reading
    config_path <- file.path(embedding_dir, 'config.json')
    jsonlite::write_json(config, config_path, pretty = TRUE, auto_unbox = TRUE)

    # save image paths and Y values
    train_csv <- paste0(embedding_dir, '/x_obs_for_train.csv')
    write.csv(X_obs, train_csv, row.names = FALSE)

    # generate X_obs
    if (embedding_obs == 'scratch') {
      # For 'scratch': train from scratch first, then use trained model
      cat("🏗️  Training model from scratch...\n")
        # Run training pipeline
      train_script <- "inst/learn_embedding/run_train_test_pipeline.py"
      args_vec <- c("--input_csv", normalizePath(train_csv),
                    "--id_col", "id",
                    "--output_dir", normalizePath(embedding_dir),
                    "--script_dir", normalizePath(script_dir),
                    "--model_name", "resnet18",
                    "--task", "regression",
                    "--epochs", "100",
                    "--batch_size", "4",
                    "--test_size", "0.5",
                    "--val_frac", "0.2",
                    "--amp",
                    "--lr", "3e-3",
                    "--use_tensorboard")

      # Use the python from the active env (auto-detected via CONDA_PREFIX)
      res <- run_python_safe(train_script, args_vec)
      cat(res, sep = "\n")  
      
    } else if (embedding_obs == 'medicalnet_ft'){
      # For 'scratch': train from scratch first, then use trained model
      cat("🏗️  Fine-tuning pretrained MedicalNet weights...\n")
        # Run training pipeline
      train_script <- "inst/learn_embedding/run_train_test_pipeline.py"
      args_vec <- c("--input_csv", normalizePath(train_csv),
                    "--id_col", "id",
                    "--output_dir", normalizePath(embedding_dir),
                    "--script_dir", normalizePath(script_dir),
                    "--model_name", "resnet18",
                    "--task", "regression",
                    "--epochs", "100",
                    "--batch_size", "4",
                    "--test_size", "0.5",
                    "--val_frac", "0.2",
                    "--amp",
                    "--lr", "3e-3",
                    "--use_tensorboard",
                    "--pretrained",
                    "--simple_head")

      # Use the python from the active env (auto-detected via CONDA_PREFIX)
      res <- run_python_safe(train_script, args_vec)
      cat(res, sep = "\n") 
    } else if(embedding_obs == 'medicalnet_ft_frozen'){
      # For 'scratch': train from scratch first, then use trained model
      cat("🏗️  Fine-tuning pretrained MedicalNet weights with frozen backbone...\n")
         # Run training pipeline
      train_script <- "inst/learn_embedding/run_train_test_pipeline.py"
      args_vec <- c("--input_csv", normalizePath(train_csv),
                    "--id_col", "id",
                    "--output_dir", normalizePath(embedding_dir),
                    "--script_dir", normalizePath(script_dir),
                    "--model_name", "resnet18",
                    "--task", "regression",
                    "--epochs", "100",
                    "--batch_size", "4",
                    "--test_size", "0.5",
                    "--val_frac", "0.2",
                    "--amp",
                    "--lr", "3e-3",
                    "--use_tensorboard",
                    "--pretrained",
                    "--simple_head",
                    "--freeze_backbone")

      # Use the python from the active env (auto-detected via CONDA_PREFIX)
      res <- run_python_safe(train_script, args_vec)
      cat(res, sep = "\n") 
    }else{
      # error
      stop("Error: embedding_obs not supported")
    }
  }

  # Load embeddings and convert to proper format
embeddings <- as.matrix(arrow::read_parquet(embeddings_path))
  index_df <- read.csv(index_path, stringsAsFactors = FALSE)

  # Convert to data frame, add IDs, and align to X_orig$id
  embedding_df <- as.data.frame(embeddings)
  X_obs <- data.frame(id = index_df[['subject_id']], embedding_df)
  # Create numeric ids for X_obs
  X_obs$id <- as.integer(X_obs$id) 
  
  #reorder X_orig, Z, Y_df to match X_obs$id (since ids now potentially only from test set of network fit)
  id_order <- match(X_obs$id, X_orig$id)
  X_orig <- X_orig[id_order, ]
  Z <- Z[id_order, , drop=FALSE]
  Y <- Y[id_order, , drop=FALSE]
  # (optional) sanity check
  stopifnot(identical(X_obs$id, X_orig$id))
}



### now with data_gen function
args <- c("/No_CI/", "1", "0", "1", "0", "fastsurfer", "scratch", "ukb_z4", "squared", "RCOT", "1")
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
beta2s=list(0.5)
idx_beta2=1
idx_sample = 1
seed=1
n_sample = list(550)

## test data generation
set.seed(seed)

##### Test DNCIT-specific embedding maps
## prepare data
path_to_ukb_data <- Sys.getenv("UKB_PATH", unset = NA)

XYZ_list <- dncitPaper::data_gen(seed=seed, idx_sample=idx_sample, n_sample=n_sample, idx_beta2=idx_beta2, 
                      beta2s=beta2s, post_non_lin=post_non_lin, eps_sigmaX=eps_sigmaX, eps_sigmaY=eps_sigmaY, 
                      eps_sigmaZ=eps_sigmaZ, embedding_orig=embedding_orig, embedding_obs=embedding_obs, 
                      confounder=confounder, g_z=g_z)

X <- as.matrix(XYZ_list[[1]])
Y <- as.matrix(XYZ_list[[2]])
Z <- as.matrix(XYZ_list[[3]])


## test for several eps_sigmaY
# for 3 seeds
idx_samples <- 1:3
n_sample <- list(460, 825, 1100)
xz_modes <- c('independent')
seeds <- 1:10
eps_sigmaY_list <- c(0, 0.05, 0.1, 0.5, 1)
for(eps_sigmaY in eps_sigmaY_list){
  for(seed in seeds){
    for(xz_mode in xz_modes){
      for(idx_sample in idx_samples){
        XYZ_list <- dncitPaper::data_gen(seed=seed, idx_sample=idx_sample, n_sample=n_sample, idx_beta2=idx_beta2, 
                            beta2s=beta2s, post_non_lin=post_non_lin, eps_sigmaX=eps_sigmaX, eps_sigmaY=eps_sigmaY, 
                            eps_sigmaZ=eps_sigmaZ, embedding_orig=embedding_orig, embedding_obs=embedding_obs, 
                            confounder=confounder, g_z=g_z, xz_mode=xz_mode)
        X <- as.matrix(XYZ_list[[1]])
        Y <- as.matrix(XYZ_list[[2]])
        Z <- as.matrix(XYZ_list[[3]])
      }
    }
  }
}