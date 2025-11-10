# load T1 paths for diagnostic
library(dplyr)
library(devtools)
load_all()
#library(reticulate)
#np <- import("numpy")
library(arrow)

# load UKB image_path.csv 
readRenviron(".Renviron")
Sys.getenv(c("UKB_PATH","IXI_PATH"))
ukb_path <- Sys.getenv("UKB_PATH", unset = NA)
t1_paths_diag <- read.csv(paste0(ukb_path, "/t1_paths_diag.csv"))

### Load embeddings
# load fastsurfer, condVAE, freesurfer embeddings
fastsurfer_embedding <- data.table::fread(paste0(ukb_path, "/ukb_fastsurfer.csv"), header=TRUE, nThread = 1)
colnames(fastsurfer_embedding)[1] <- "id"
condVAE_embedding <- data.table::fread(paste0(ukb_path, "/ukb_condVAE.csv"), header=FALSE, nThread = 1)
colnames(condVAE_embedding)[1] <- "id"
freesurfer_embedding <- data.table::fread(paste0(ukb_path, "/ukb_freesurfer.csv"), header=TRUE, nThread = 1)
colnames(freesurfer_embedding)[1] <- "id"

# Extract pretrained MedicalNet embeddings
out_dir <- paste0(ukb_path, "/medicalnet_embeddings")
#system2("python", c("inst/learn_embedding/extract_embeddings.py",
#                    "--pretrained", "--model_name", "resnet18",
#                    "--input_csv", paste0(ukb_path, "/t1_paths.csv"),
#                    "--output_dir", out_dir, 
#                    "--batch_size", "16",
#                    "--num_workers", "4"))
medicalnet_emb <- as.matrix(arrow::read_parquet(paste0(out_dir, "/embeddings.parquet")))
medicalnet_idx <- fread(paste0(out_dir, "/embeddings_index.csv"))
medicalnet_embedding <- cbind(id = medicalnet_idx$subject_id, as.data.frame(medicalnet_emb))

# Find common IDs across all datasets
# Assuming Y is stored in a dataframe/data.table with an 'id' column
common_ids <- Reduce(intersect, list(
  t1_paths_diag$id,
  fastsurfer_embedding$id,  # or $id depending on column name
  condVAE_embedding$id,
  freesurfer_embedding$id,
  medicalnet_embedding$id
))
cat("Total common IDs across all datasets:", length(common_ids), "\n")
common_ids <- sort(common_ids)

# Filter and reorder each dataset
t1_paths_diag <- t1_paths_diag[match(common_ids, t1_paths_diag$id), ]
fastsurfer_emb <- fastsurfer_embedding[match(common_ids, fastsurfer_embedding$id), ]  # adjust column name
condVAE_emb <- condVAE_embedding[match(common_ids, condVAE_embedding$id), ]
freesurfer_emb <- freesurfer_embedding[match(common_ids, freesurfer_embedding$id), ]
medicalnet_emb <- medicalnet_embedding[match(common_ids, medicalnet_embedding$id), ]

write.csv(t1_paths_diag, file.path(ukb_path, 't1_paths_diag_cvae.csv'), row.names = FALSE)

# Extract embeddings from trained checkpoint
current_y_gen_dir <- paste0(ukb_path, "/No_CI/550/1")
checkpoint_path <- paste0(current_y_gen_dir, "/scratch/best_model.ckpt")
out_dir_diagnostic <- paste0(current_y_gen_dir, "/diagnostic")
dir.create(out_dir_diagnostic, recursive = TRUE)
system2("python", c("inst/learn_embedding/extract_embeddings.py",
                    "--checkpoint", checkpoint_path,
                    "--input_csv", paste0(ukb_path, "/t1_paths_diag_cvae.csv"),
                    "--output_dir", out_dir_diagnostic,
                    "--config_name", "scratch", 
                    "--batch_size", "16",
                    "--num_workers", "4",
                    "--amp"))
trained_emb <- as.matrix(arrow::read_parquet(paste0(out_dir_diagnostic, "/scratch/embeddings.parquet")))
trained_idx <- fread(paste0(out_dir_diagnostic, "/scratch/embeddings_index.csv"))
trained_embedding <- cbind(id = trained_idx$subject_id, as.data.frame(trained_emb))




### gen Y for current setting
# set up as in sim_run_settings.R and sim_ukb_brainmri.R
args <- c("/No_CI/", "1", "0", "0.05", "0", "fastsurfer", "scratch", "ukb_z4", "squared", "RCOT", "1")
n_cits <- 1
cit <- c(args[10])
print(args)
post_non_lin=as.numeric(args[2])
eps_sigmaX=as.numeric(args[3])
eps_sigmaY=as.numeric(args[4])
eps_sigmaZ=as.numeric(args[5])
embedding_orig=args[6]
embedding_obs=args[7]
confounder=args[8]                                        
g_z=args[9]
if(args[1] == "/No_CI/"){
  beta2s=list(1)
  idx_beta2=1
}
seed = 1
set.seed(seed)
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
is_binary <- function(x) {
  unique_values <- unique(x)
  length(unique_values) == 2
}

Z <- load_Z(ukb_path,confounder)

epsZ <- matrix(stats::rnorm((nrow(Z)*ncol(Z)), 0, eps_sigmaZ), nrow=nrow(Z), ncol=ncol(Z))
Z <- as.data.frame(Z)+epsZ
Z <- Z[match(common_ids, trained_embedding$id), ]

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
  site_sum <- sum(as.data.frame(Z)[site_columns])
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
# scale fastsurfer embeddings (X_orig)
fastsurfer_emb <- as.data.frame(fastsurfer_emb)
fastsurfer_emb[,c(-1)] <- scale(fastsurfer_emb[,c(-1)])

if(is.null(beta2s)){
  Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z)
}else{
  Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, X=as.matrix(fastsurfer_emb[,c(-1)]), beta2s=beta2s, idx_beta2=idx_beta2, post_non_lin=post_non_lin, g_z=g_z)
}

Y_id <- data.frame(id = fastsurfer_emb$id, Y = Y[,1])
write.csv(Y_id, file.path(out_dir_diagnostic, 'Y_diag.csv'), row.names = FALSE)

#### predict Y using glmnet
library(glmnet)

ridge_lasso_glmnet <- function(
  X, y,
  alpha = 1,                 # 1 = lasso, 0 = ridge, (0,1) = elastic net
  test_prop = 0.2,
  nfolds = 10,
  lambda_choice = c("min", "1se"),
  seed = 42,
  standardize = TRUE,
  intercept = TRUE,
  return_model = FALSE
) {
  stopifnot(is.numeric(y))
  X <- as.matrix(X)
  stopifnot(nrow(X) == length(y))

  set.seed(seed)
  n <- nrow(X)
  idx_test  <- sample.int(n, size = max(1, round(test_prop * n)))
  idx_train <- setdiff(seq_len(n), idx_test)

  Xtr <- X[idx_train, , drop = FALSE]
  ytr <- y[idx_train]
  Xte <- X[idx_test,  , drop = FALSE]
  yte <- y[idx_test]

  # cross-validated glmnet
  suppressPackageStartupMessages(require(glmnet))
  cvfit <- cv.glmnet(
    x = Xtr, y = ytr,
    family = "gaussian",
    alpha = alpha,
    nfolds = nfolds,
    standardize = standardize,
    intercept = intercept
  )

  lambda_choice <- match.arg(lambda_choice)
  s_use <- if (lambda_choice == "min") cvfit$lambda.min else cvfit$lambda.1se

  pred <- drop(predict(cvfit, newx = Xte, s = s_use))
  mse  <- mean((yte - pred)^2)
  r2   <- 1 - sum((yte - pred)^2) / sum((yte - mean(yte))^2)

  out <- list(
    mse_test = as.numeric(mse),
    r2_test  = as.numeric(r2),
    lambda   = s_use,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    alpha = alpha,
    idx_train = idx_train,
    idx_test = idx_test,
    y_test = yte,
    yhat_test = pred
  )
  if (return_model) out$cvfit <- cvfit
  out
}

fastsurfer_results <- ridge_lasso_glmnet(X=fastsurfer_emb, y=Y, lambda_choice = "1se")
condVAE_results <- ridge_lasso_glmnet(X=condVAE_emb, y=Y, lambda_choice = "1se")
freesurfer_results <- ridge_lasso_glmnet(freesurfer_emb, Y, lambda_choice = "1se")
medicalnet_results <- ridge_lasso_glmnet(medicalnet_emb, Y, lambda_choice = "1se")
trained_results <- ridge_lasso_glmnet(trained_emb, Y, lambda_choice = "1se")

print(paste0("FastSurfer: R^2 = ", fastsurfer_results$r2_test, " MSE = ", fastsurfer_results$mse_test))
print(paste0("condVAE: R^2 = ", condVAE_results$r2_test, " MSE = ", condVAE_results$mse_test))
print(paste0("Freesurfer: R^2 = ", freesurfer_results$r2_test, " MSE = ", freesurfer_results$mse_test))
print(paste0("MedicalNet: R^2 = ", medicalnet_results$r2_test, " MSE = ", medicalnet_results$mse_test))
print(paste0("Trained: R^2 = ", trained_results$r2_test, " MSE = ", trained_results$mse_test))



#### via auto_diagnostic function
idx_samples <- 1:4
n_sample = list(460, 1100, 5000, 10000)
xz_modes <- c('independent')
seeds <- c(301, 303, 312:314, 323:325)
eps_sigmaY_list <- c(1)
embedding_obs <- c('medicalnet_ft_frozen', 'medicalnet_ft', 'scratch')
Y_age <- FALSE
standardize_ridge_lasso <- TRUE
for(eps_sigmaY in eps_sigmaY_list){
  for(seed in seeds){
    # Cache baseline results per seed (constant across sample sizes)
    baseline_cache <- NULL
    for(xz_mode in xz_modes){
      for(idx_sample in idx_samples){
        if(idx_sample == 1){
          # Compute baseline results for first sample size and cache them
          results <- auto_diagnostic(
            experiment_dir = paste0("/sc/home/marco.simnacher/ukbiobank/data/No_CI/", n_sample[idx_sample], "/", seed, "/eps_sigmaY=", as.character(eps_sigmaY)),
            embedding_obs = embedding_obs,
            seed = seed,
            debug_Y = FALSE,
            lambda_choice = "min",
            Y_age = Y_age,
            standardize_ridge_lasso = standardize_ridge_lasso
          )
          baseline_cache <- results$baseline_results
        } else {
          # Reuse cached baseline results for subsequent sample sizes
          results <- auto_diagnostic(
            experiment_dir = paste0("/sc/home/marco.simnacher/ukbiobank/data/No_CI/", n_sample[idx_sample], "/", seed, "/eps_sigmaY=", as.character(eps_sigmaY)),
            embedding_obs = embedding_obs,
            seed = seed,
            debug_Y = FALSE,
            baseline_results_cached = baseline_cache,
            lambda_choice = "min",
            Y_age = Y_age,
            standardize_ridge_lasso = standardize_ridge_lasso
          )
        }
      }
    }
  }
}
