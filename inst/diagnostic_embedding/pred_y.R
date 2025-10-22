# load T1 paths for diagnostic
library(data.table)
library(dplyr)
library(devtools)
load_all()
library(reticulate)
np <- import("numpy")

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
system2("python", c("inst/learn_embedding/extract_embeddings.py",
                    "--pretrained", "--model_name", "resnet18",
                    "--input_csv", paste0(ukb_path, "/t1_paths_diag.csv"),
                    "--output_dir", out_dir, 
                    "--batch_size", "16",
                    "--num_workers", "4"))
medicalnet_emb <- np$load(paste0(out_dir, "/embeddings.npy"))
medicalnet_idx <- fread(paste0(out_dir, "/embeddings_index.csv"))
medicalnet_embedding <- cbind(id = medicalnet_idx$subject_id, as.data.frame(medicalnet_emb))

# Extract embeddings from trained checkpoint
checkpoint_path <- paste0(ukb_path, "/t1/embeddings/550_scratch_ukb_z4_linear/1/best_model.ckpt")
out_dir_trained <- paste0(ukb_path, "/trained_model_embeddings")
system2("python", c("inst/learn_embedding/extract_embeddings.py",
                    "--checkpoint", checkpoint_path,
                    "--input_csv", paste0(ukb_path, "/t1_paths_diag.csv"),
                    "--output_dir", out_dir_trained,
                    "--config_name", "trained", 
                    "--batch_size", "16",
                    "--num_workers", "4"))
trained_emb <- np$load(paste0(out_dir_trained, "/trained/embeddings.npy"))
trained_idx <- fread(paste0(out_dir_trained, "/trained/embeddings_index.csv"))
trained_embedding <- cbind(id = trained_idx$subject_id, as.data.frame(trained_emb))

Y <- read.csv(paste0(ukb_path, "/ukb_z1_age.csv"))

# Find common IDs across all datasets
# Assuming Y is stored in a dataframe/data.table with an 'id' column
common_ids <- Reduce(intersect, list(
  t1_paths_diag$id,
  fastsurfer_embedding$id,  # or $id depending on column name
  #condVAE_embedding$id,
  freesurfer_embedding$id,
  medicalnet_embedding$id,
  trained_embedding$id,
  Y$id  # adjust column name as needed
))

cat("Total common IDs across all datasets:", length(common_ids), "\n")

# Filter and sort all datasets by common IDs
# Sort common_ids first for consistent ordering
common_ids <- sort(common_ids)

# Filter and reorder each dataset
t1_paths_diag <- t1_paths_diag[match(common_ids, t1_paths_diag$id), ]
fastsurfer_emb <- fastsurfer_embedding[match(common_ids, fastsurfer_embedding$id), ]  # adjust column name
condVAE_emb <- condVAE_embedding[match(common_ids, condVAE_embedding$id), ]
freesurfer_emb <- freesurfer_embedding[match(common_ids, freesurfer_embedding$id), ]
medicalnet_emb <- medicalnet_embedding[match(common_ids, medicalnet_embedding$id), ]
trained_emb <- trained_embedding[match(common_ids, trained_embedding$id), ]
Y_aligned <- Y[match(common_ids, Y$id), ]

# Verify alignment (all should return TRUE)
cat("Checking alignment:\n")
cat("  t1_paths_diag IDs match:", all(t1_paths_diag$id == common_ids), "\n")
cat("  fastsurfer IDs match:", all(fastsurfer_emb$eid == common_ids), "\n")
#cat("  condVAE IDs match:", all(condVAE_emb$id == common_ids), "\n")
cat("  freesurfer IDs match:", all(freesurfer_emb$id == common_ids), "\n")
cat("  medicalnet IDs match:", all(medicalnet_emb$id == common_ids), "\n")
cat("  trained IDs match:", all(trained_emb$id == common_ids), "\n")
cat("  Y IDs match:", all(Y_aligned$id == common_ids), "\n")

# ---- train / test split ----
n  <- length(common_ids)
idx_test  <- sample.int(n, size = round(0.2*n))   # 20% hold-out
idx_train <- setdiff(seq_len(n), idx_test)

fastsurfer_tr <- fastsurfer_emb[idx_train, -1, drop = FALSE]
fastsurfer_te <- fastsurfer_emb[idx_test, -1, drop = FALSE]
freesurfer_tr <- freesurfer_emb[idx_train, -1, drop = FALSE]
freesurfer_te <- freesurfer_emb[idx_test, -1, drop = FALSE]
medicalnet_tr <- medicalnet_emb[idx_train, -1, drop = FALSE]
medicalnet_te <- medicalnet_emb[idx_test, -1, drop = FALSE]
trained_tr <- trained_emb[idx_train, -1, drop = FALSE]
trained_te <- trained_emb[idx_test,  -1, drop = FALSE]
ytr <- Y_aligned[idx_train,2]
yte <- Y_aligned[idx_test,2]

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

fastsurfer_results <- ridge_lasso_glmnet(X=fastsurfer_emb, y=Y_aligned$age, lambda_choice = "1se")
freesurfer_results <- ridge_lasso_glmnet(freesurfer_emb, Y_aligned$age, lambda_choice = "1se")
medicalnet_results <- ridge_lasso_glmnet(medicalnet_emb, Y_aligned$age, lambda_choice = "1se")
trained_results <- ridge_lasso_glmnet(trained_emb, Y_aligned$age, lambda_choice = "1se")

print(paste0("FastSurfer: ", fastsurfer_results$r2_test))
print(paste0("Freesurfer: ", freesurfer_results$r2_test))
print(paste0("MedicalNet: ", medicalnet_results$r2_test))
print(paste0("Trained: ", trained_results$r2_test))