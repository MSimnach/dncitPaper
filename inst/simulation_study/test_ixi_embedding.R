##### Embedding maps trained on IXI dataset
# Helper function for safe Python execution
run_python_safe <- function(script_path, args = "") {
  helper_script <- "/home/RDC/simnacma/Coding/dncitPaper/inst/simulation_study/run_python_with_conda_libs.sh"
  cmd <- sprintf("%s %s %s", helper_script, script_path, args)
  system(cmd, intern = TRUE)
}

setwd("./dncitPaper")
library(devtools)
load_all()
# Load embeddings using reticulate (numpy only - should work fine)
library(reticulate)
np <- import("numpy")

# load IXI image_path.csv after running prepare_ixi_from_tar_xls.py
readRenviron("/home/RDC/simnacma/Coding/dncitPaper/.Renviron")
Sys.getenv(c("UKB_PATH","IXI_PATH"))
ixi_path <- Sys.getenv("IXI_PATH", unset = NA)
script_dir <- "/home/RDC/simnacma/Coding/dncitPaper/inst/learn_embedding"

ixi_data <- read.csv(paste0(ixi_path, "/t1/csv/all_labeled.csv"))

# Example brain MRI plot
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

# Audit results of IXI data (raw images, for orientation, scaling, etc.)
ixi_audit <- read.csv(paste0(ixi_path, "/t1/qc_summary.csv"))

# y check
ixi_y <- read.csv(paste0(ixi_path, "/t1/csv/all_labeled.csv"))
ixi_y$y <10

# regression results check
ixi_regression_check <- read.csv(paste0(ixi_path, "/t1/results/analysis/per_observation_results.csv"))
ixi_regression_check[order(ixi_regression_check$squared_error),]
ixi_y[ixi_y$y>80,]
ixi_y[ixi_y$id=="IXI045",]


# paths
res_dir <- "/home/RDC/simnacma/H:/simnacma/CITs/Application/ixi/t1/results/scratch"
pred_path  <- file.path(res_dir, "test_predictions.npy")
index_path <- file.path(res_dir, "test_embeddings_index.csv")
test_csv   <- "/home/RDC/simnacma/H:/simnacma/CITs/Application/ixi/t1/csv/test_labeled.csv"

# load predictions (ensure 1-D numeric vector) - using same np from above
ixi_test_preds <- as.numeric(np$load(pred_path))
length(ixi_test_preds)

# load index and labels
ixi_test_index  <- read.csv(index_path, stringsAsFactors = FALSE)
ixi_test_labels <- read.csv(test_csv, stringsAsFactors = FALSE)

names(ixi_test_index)  # inspect columns
names(ixi_test_labels) # typically c("id","path","y")

# The index CSV written by our scripts usually has columns: "row","subject_id"
# Align by row order and then attach the subject id:
stopifnot(nrow(ixi_test_index) == length(ixi_test_preds))
df_pred <- data.frame(
  row  = ixi_test_index$row,
  id   = if ("subject_id" %in% names(ixi_test_index)) ixi_test_index$subject_id else ixi_test_index$id,
  yhat = ixi_test_preds,
  stringsAsFactors = FALSE
)

# Now merge with the labeled test CSV by 'id'
key_col <- if ("id" %in% names(ixi_test_labels)) "id" else "subject_id"
merged <- merge(df_pred, ixi_test_labels[, c(key_col, "y")],
                by.x = "id", by.y = key_col, all.x = TRUE)

head(merged)


############## set up as in sim_run_settings.R and sim_ukb_brainmri.R
args <- c("/CI/", "1", "0", "1", "0", "medicalnet", "scratch", "ukb_z4", "linear", "RCOT", "1")
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
beta2s=NULL
idx_beta2=NULL
idx_sample = 1
seed=1
n_sample = list(550)
## test data generation
set.seed(seed)

##### Test DNCIT-specific embedding maps
## prepare data
path_to_ixi_data <- Sys.getenv("IXI_PATH", unset = NA)

# data loading functions
load_X_orig <- function(path_to_ixi_data, embedding_orig){
  if(embedding_orig == 'medicalnet'){
    X <- data.table::fread(paste0(path_to_ixi_data, '/t1/csv/all_labeled.csv'), header=TRUE, nThread = 1)
  }
  return(X)
}

load_X_obs <- function(path_to_ixi_data,embedding_obs, embedding_orig, X_orig, eps_sigmaX){
  if(embedding_obs %in% c('medicalnet', 'scratch', 'ft_head_only', 'ft_full')){
    X_obs <- data.table::fread(paste0(path_to_ixi_data, '/t1/csv/all_labeled.csv'), header=TRUE, nThread = 1)
  }
  return(X_obs)
}

load_Z <- function(path_to_ixi_data,confounder){
  if(confounder=='ukb_z4'){
    Z <- data.table::fread(paste0(path_to_ixi_data, "/ixi_confounder.csv"), header=TRUE, nThread = 1)
  }
}

Z <- load_Z(path_to_ixi_data,confounder)
X_orig <- load_X_orig(path_to_ixi_data,embedding_orig)
X_obs <- load_X_obs(path_to_ixi_data,embedding_obs, embedding_orig, X_orig, eps_sigmaX)

# take only common IDs of X, Z
colnames(X_orig)[1] <- 'id'
colnames(X_obs)[1] <- 'id'
colnames(Z)[1] <- 'id'
numeric_ids <- as.integer(gsub("IXI", "", X_orig$id))
X_orig$id <- numeric_ids
numeric_ids <- as.integer(gsub("IXI", "", X_obs$id))
X_obs$id <- numeric_ids
X_orig <- as.data.frame(X_orig)
X_obs <- as.data.frame(X_obs)
Z <- as.data.frame(Z)

Z <- Z[!duplicated(Z$id), ]
common_ids <- Reduce(intersect, list(X_orig$id, X_obs$id, Z$id))

subset_X_obs <- X_obs[X_obs$id %in% common_ids, ]
subset_X_orig <-  X_orig[X_orig$id %in% common_ids, ]
subset_Z <- Z[Z$id %in% common_ids, ]
#check if ids are equal
stopifnot(all.equal(subset_X_orig[,1], subset_Z[,1]))

#sample rows
idx <- sample(1:nrow(subset_Z), n_sample[[idx_sample]])
X_orig <- subset_X_orig[idx,]
X_obs <- subset_X_obs[idx,]
Z <- subset_Z[idx,]

epsZ <- matrix(stats::rnorm((nrow(Z)*ncol(Z)), 0, eps_sigmaZ), nrow=nrow(Z), ncol=ncol(Z))
Z <- Z+epsZ

#Standardize
# scale only continuous confounders
Z[,c(-1)] <- Z[,c(-1)] %>% dplyr::mutate(dplyr::across(dplyr::where(function(x) !is_binary(x)), scale))

if(embedding_orig == 'medicalnet'){
  #get embeddings of X_orig and order w.r.t. idx, i.e. Z
  medicalnet_dir <- paste0(path_to_ixi_data, '/t1/embeddings/', embedding_orig)
  embeddings_orig <- np$load(file.path(medicalnet_dir, 'embeddings.npy'))
  embedding_orig_index  <- read.csv(file.path(medicalnet_dir, 'embeddings_index.csv'), stringsAsFactors = FALSE)
  # use only embeddings with index in X_orig$id (remove IXI from id in embedding_orig_index$subject_id and remove 0's if number starts with 0)
  embedding_orig_index$subject_id_num <- as.integer(sub("^IXI", "", embedding_orig_index$subject_id, ignore.case = TRUE))
  # Row indices in embedding table that match X_orig$id (in X_orig's order)
  order_idx <- match(X_orig$id, embedding_orig_index$subject_id_num)
  # Reorder embeddings to X_orig$id order
  embeddings_orig_aligned <- embeddings_orig[order_idx, , drop = FALSE]

  X_orig <- data.frame(id = embedding_orig_index$subject_id_num[order_idx], embeddings_orig_aligned)
}

# scale X_orig
X_orig[,c(-1)] <- scale(X_orig[,c(-1)])

# generate Y
if(is.null(beta2s)){
  Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z)
}else{
  Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, X=X_orig[,c(-1)], beta2s=beta2s, idx_beta2=idx_beta2, post_non_lin=post_non_lin, g_z=g_z)
}


if(embedding_obs %in% c('fastsurfer', 'condVAE', 'latentDiffusion', 'freesurfer')){
  X_obs[,c(-1)] <- scale(X_obs[,c(-1)])
} else{
  # Reorder X_obs to match X_orig,Z,Y and assign Y to X_obs
  X_obs <- X_obs[match(X_orig$id, X_obs$id), ]
  X_obs$y <- Y[,1]

  # Setup paths for train results and embeddings
  input_csv <- paste0(path_to_ixi_data, '/t1/csv/all_labeled.csv')
  embedding_dir <- paste0(path_to_ixi_data, '/t1/embeddings/', embedding_obs, '/', seed)
  train_output_dir <- paste0(path_to_ixi_data, '/t1/results/', embedding_obs)

  embeddings_path <- file.path(embedding_dir, 'test_embeddings.npy')
  index_path <- file.path(embedding_dir, 'test_embeddings_index.csv')
  preds_path <- file.path(embedding_dir, 'test_predictions.npy')

  ## Train or skip if embeddings already exist
  if (dir.exists(embedding_dir)){
    cat("✅ Embedding directory already exists. Skipping training...\n")
  } else {
    dir.create(embedding_dir, recursive = TRUE)
    # save Y and data gen config
    # Save configuration
    config <- list(
      args = args,
      config_names = c("path", "post_non_lin", "eps_sigmaX", "eps_sigmaY", "eps_sigmaZ", 
                      "embedding_orig", "embedding_obs", "confounder", "g_z", "cit"),
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
    # Save Y values
    Y_path <- file.path(embedding_dir, 'Y.csv')
    Y_df <- data.frame(
      id = X_orig$id,
      Y = Y[,1]
    )
    write.csv(Y_df, Y_path, row.names = FALSE)

    train_csv <- paste0(embedding_dir, '/x_obs_for_train.csv')
    write.csv(X_obs, train_csv, row.names = FALSE)
    
    # generate X_obs
    if (embedding_obs == 'scratch') {
      # For 'scratch': train from scratch first, then use trained model
      cat("🏗️  Training model from scratch...\n")
        # Run training pipeline
      train_script <- "/home/RDC/simnacma/Coding/dncitPaper/inst/learn_embedding/run_train_test_pipeline.py"
      train_args <- sprintf("--input_csv %s --id_col id --output_dir %s --script_dir %s --model_name resnet18 --task regression --epochs 30 --batch_size 16 --test_size 0.5 --val_frac 0.2 --amp --num_workers 4 --lr 8e-4",
                          train_csv, train_output_dir, script_dir)
      
      cat("Training with args:", train_args, "\n")
      result <- run_python_safe(train_script, train_args)
      
    } else if (embedding_obs == 'medicalnet_ft'){
      # For 'medicalnet_ft': fine-tune pretrained medicalnet weights, then use fine-tuned model
      cat("🏗️  Fine-tune model with pretrained medicalnet weights...\n")
      # Run training pipeline
      train_script <- "/home/RDC/simnacma/Coding/dncitPaper/inst/learn_embedding/run_train_test_pipeline.py"
      train_args <- sprintf("--input_csv %s --id_col id --output_dir %s --script_dir %s --model_name resnet18 --task regression --epochs 30 --batch_size 16 --test_size 0.5 --val_frac 0.2 --amp --num_workers 4 --lr 8e-4 --pretrained --simple_head",
                          train_csv, train_output_dir, script_dir)
      
      cat("Training with args:", train_args, "\n")
      result <- run_python_safe(train_script, train_args)
    } else{
      # For pretrained models: use extraction script
      embedding_dir <- paste0(path_to_ixi_data, '/t1/embeddings/', embedding_obs)
      embeddings_path <- file.path(embedding_dir, 'embeddings.npy')
      index_path <- file.path(embedding_dir, 'embeddings_index.csv')
      
      # Create embedding directory
      dir.create(embedding_dir, recursive = TRUE)
      
      # Determine model name
      model_name <- switch(embedding_obs,
                          'medicalnet' = 'resnet18',
                          'ft_head_only' = 'resnet18',
                          'ft_full' = 'resnet18',
                          'resnet18')
      
      # Run extraction using helper script
      python_script <- "/home/RDC/simnacma/Coding/dncitPaper/inst/learn_embedding/run_ixi_extraction.py"
      batch_size <- 16
      args <- sprintf("--input_csv %s --output_dir %s --model_name %s --batch_size %d",
                      input_csv, embedding_dir, model_name, batch_size)
      
      cat("Running:", model_name, "on", input_csv, "\n")
      result <- run_python_safe(python_script, args)
    }
    # copy results from train_output_dir to embedding_dir
    file.copy(from = file.path(train_output_dir, 'test_embeddings.npy'), to = embeddings_path)
    file.copy(from = file.path(train_output_dir, 'test_embeddings_index.csv'), to = index_path)
    file.copy(from = file.path(train_output_dir, 'test_predictions.npy'), to = preds_path)
  }

  # Load embeddings and convert to proper format
  embeddings <- np$load(embeddings_path)
  index_df <- read.csv(index_path, stringsAsFactors = FALSE)

  # Convert to data frame, add IDs, and align to X_orig$id
  embedding_df <- as.data.frame(embeddings)
  X_obs <- data.frame(id = index_df[['subject_id']], embedding_df)
  # Create numeric ids for X_obs
  X_obs$id <- as.integer(sub("^IXI", "", X_obs$id, ignore.case = TRUE))
  
  #reorder X_orig, Z, Y_df to match X_obs$id (since ids now potentially only from test set of network fit)
  id_order <- match(X_obs$id, X_orig$id)
  X_orig <- X_orig[id_order, ]
  Z <- Z[id_order, , drop=FALSE]
  Y <- Y[id_order, , drop=FALSE]
  # (optional) sanity check
  stopifnot(identical(X_obs$id, X_orig$id))
}

row.names(Z) <- 1:nrow(Z)
row.names(X_obs) <- 1:nrow(X_obs)
row.names(X_orig) <- 1:nrow(X_orig)
row.names(Y) <- 1:nrow(Y)

XYZ_list <- list(X_obs[,c(-1)],Y,Z[,c(-1)])



#### CIT
X <- as.matrix(XYZ_list[[1]])
Y <- as.matrix(XYZ_list[[2]])
Z <- as.matrix(XYZ_list[[3]])

if (args[10] == 'RCOT'){
                                                       cit_params <- list(cit='RCOT', params_cit=list(seed=as.numeric(args[11])))
                                                     }else if(args[10] == 'CMIknn'){
                                                       cit_params <- list(cit='cmiknn', params_cit=list())
                                                     }else if(args[10] == 'kpc_graph'){
                                                      if (args[11]=='1'){
                                                          k = kernlab::vanilladot()
                                                      }else if (args[11]=='2') {
                                                          k = kernlab::rbfdot(1/(2*stats::median(stats::dist(X))^2))
                                                      }else if (args[11]=='3') {
                                                          k = kernlab::laplacedot(1/(2*stats::median(stats::dist(X))^2))
                                                      }else if (args[11]=='4') {
                                                          k = kernlab::tanhdot()
                                                      }
                                                      if(args[8]=='ukb_z1'){
                                                        model_formula_YZ <- "V1~1+s(V2)"
                                                      }else if(args[8]=='ukb_z2'){
                                                        model_formula_YZ <- 'V1~1+s(V2)+s(V3)'
                                                      }else if(args[8]=='ukb_z4'){
                                                        n_covs <- ncol(Z)
                                                        lin_covs <- paste0("V", seq(4, n_covs+1))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste0('V1~1+s(V2)+s(V3)+', lin_covs_string)
                                                      }else if(args[8]=='ukb_z6'){
                                                        n_covs <- ncol(Z)
                                                        date_diff <- as.character(n_covs)
                                                        qc <- as.character(n_covs+1)
                                                        lin_covs <- paste0("V", seq(4, n_covs-1))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+', lin_covs_string, sep="")
                                                      }else if(args[8]=='ukb_z10'){
                                                        n_covs <- ncol(Z)
                                                        date_diff <- as.character(n_covs-4)
                                                        qc <- as.character(n_covs-3)
                                                        head_loc_1 <- as.character(n_covs-2)
                                                        head_loc_2 <- as.character(n_covs-1)
                                                        head_loc_3 <- as.character(n_covs)
                                                        head_loc_4 <- as.character(n_covs+1)
                                                        lin_covs <- paste0("V", seq(4, n_covs-5))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+s(V', head_loc_1, ', k=3)+s(V', head_loc_2, ', k=3)+s(V', head_loc_3, ', k=3)+s(V', head_loc_4, ', k=3)+', lin_covs_string, sep="")
                                                      }else if(args[8]=='ukb_z15'){
                                                        n_covs <- ncol(Z)
                                                        date_diff <- as.character(n_covs-9)
                                                        qc <- as.character(n_covs-8)
                                                        head_loc_1 <- as.character(n_covs-7)
                                                        head_loc_2 <- as.character(n_covs-6)
                                                        head_loc_3 <- as.character(n_covs-5)
                                                        head_loc_4 <- as.character(n_covs-4)
                                                        gene_1 <- as.character(n_covs-3)
                                                        gene_2 <- as.character(n_covs-2)
                                                        gene_3 <- as.character(n_covs-1)
                                                        gene_4 <- as.character(n_covs)
                                                        gene_5 <- as.character(n_covs+1)
                                                        lin_covs <- paste0("V", seq(4, n_covs-10))
                                                        lin_covs_string <- paste(lin_covs, collapse = "+")
                                                        model_formula_YZ <- paste('V1~1+s(V2)+s(V3)+s(V', date_diff, ')+s(V', qc, ')+s(V', head_loc_1, ', k=3)+s(V', head_loc_2, ', k=3)+s(V', head_loc_3, ', k=3)+s(V', head_loc_4, ', k=3)+s(V', gene_1, ', k=3)+s(V', gene_2, ', k=3)+s(V', gene_3, ', k=3)+s(V', gene_4, ', k=3)+s(V', gene_5, ', k=3)+', lin_covs_string, sep="")
                                                      }
                                                      cit_params <- list(cit='cpt_kpc', params_cit=list(k=k, Knn = as.numeric(args[12]), model.formula.YZ=model_formula_YZ))
                                                     }else if(args[10]=='FCIT'){
                                                        cit_params <- list(cit='fcit')
                                                     }else if(args[10]=='cpi'){
                                                        cit_params <- list(cit='cpi')
                                                     }else if(args[10]=='comets_pcm'){
                                                       cit_params <- list(cit='comets')
                                                     }else if(args[10]=='comets_gcm'){
                                                       cit_params <- list(cit='comets', params_cit=list(method='gcm', alternative='less'))
                                                     }else if(args[10]=='ccit'){
                                                       cit_params <- list(cit='ccit', params_cit=list(nthread=as.integer(1)))
                                                     }else if(args[10]=='pred_cit'){
                                                        min_samples <- min(unlist(n_sample))
                                                        max_samples <- max(unlist(n_sample))
                                                        current_sample <- n_sample[[idx_sample]]
                                                        term_time <- round(exp(1.5)+(current_sample-min_samples)/(max_samples-min_samples)*(exp(2.25)-exp(1.5))/3)
                                                        cit_params <- list(cit='pred_cit', params_cit=list(term_time = term_time))
                                                     }else if(args[10]=='WALD'){
                                                        cit_params <- list(cit='wald')
                                                     }

tmp <- DNCIT::DNCIT(X, Y, Z, embedding_map_with_parameters = 'feature_representations',
                                                           cit_with_parameters = cit_params)


