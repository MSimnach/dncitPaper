## Data generation
#' Title
#'
#' @param seed random seed for current data generation
#' @param idx_sample index of current sample size
#' @param n_sample list of sample sizes
#' @param idx_beta2 index of current beta2
#' @param beta2s list of beta2s
#' @param n sample size
#' @param post_non_lin post non-linearity function
#' @param eps_sigmaX noise level for X
#' @param eps_sigmaY noise level for Y
#' @param eps_sigmaZ noise level for Z
#' @param embedding_orig embedding map used to obtain original feature representations (used to simulate conditional association or CI)
#' @param embedding_obs embedding map used to obtain observed feature representations (used for CI testing)
#' @param confounder confounder used to simulate Y
#' @param g_z confounder response relationship
#'
#' @return list of X_obs, Y, Z
#' @export
data_gen <- function(seed, idx_sample=NULL, n_sample=NULL, idx_beta2=NULL, beta2s=NULL, post_non_lin=4, eps_sigmaX=0, eps_sigmaY=1, eps_sigmaZ=0, embedding_orig='fastsurfer',
                     embedding_obs='fastsurfer', confounder='AS', g_z='linear', debug=FALSE, xz_mode='Sigma=I_p'){
  path_to_ukb_data <- Sys.getenv("UKB_PATH", unset = NA)

  if (is.na(path_to_ukb_data) || path_to_ukb_data == "") {
    stop("❌ Environment variable UKB_PATH is not set. Please define it in your .Renviron.")
  }
  set.seed(seed)
  training_time <- 0

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

  subset_X_obs <- X_obs[match(common_ids, X_obs$id), ]
  subset_X_orig <-  X_orig[match(common_ids, X_orig$id), ]
  subset_Z <- Z[match(common_ids, Z$id), ]
  #check if ids are equal
  stopifnot(all.equal(subset_X_orig[,1], subset_Z[,1]))

  #sample rows
  idx <- sample(1:nrow(subset_Z), n_sample[[idx_sample]])
  X_obs <- subset_X_obs[idx,]
  X_orig <- subset_X_orig[idx,]
  Z <- subset_Z[idx, , drop=FALSE]

  epsZ <- matrix(stats::rnorm((nrow(Z)*(ncol(Z)-1)), 0, eps_sigmaZ), nrow=nrow(Z), ncol=ncol(Z)-1)
  Z[,c(-1)] <- Z[,c(-1), drop=FALSE]+epsZ

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
  if(isTRUE(site_sum == nrow(Z))){
    for(site in site_columns){
      if(is_binary(Z[[site]])){
        Z <- Z[, !names(Z) %in% site]
        break
      }
    }
  }

  #Standardize
  # scale only continuous confounders
  Z[,c(-1)] <- Z[,c(-1), drop=FALSE] %>% dplyr::mutate(dplyr::across(dplyr::where(function(x) !is_binary(x)), scale))

  X_orig[,c(-1)] <- scale(X_orig[,c(-1), drop=FALSE])
  #X_obs <- scale(X_obs)

  if(is.null(beta2s)){
    # Set seed specifically for Y generation (deterministic across different sample sizes)
    set.seed(seed + 999999)  # offset to separate from data sampling seed
    if (debug) {
      Y <- y_from_xz(Z[,c(-1), drop=FALSE], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z, debug=debug)
    } else {
      Y <- y_from_xz(Z[,c(-1), drop=FALSE], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z)
    }
    is_ci <- "CI"
  }else{
    # Set seed specifically for Y generation (deterministic across different sample sizes)
    set.seed(seed + 999999)  # offset to separate from data sampling seed
    if (debug) {
      Y <- y_from_xz(Z[,c(-1), drop=FALSE], eps_sigmaY, X=as.matrix(X_orig[,c(-1), drop=FALSE]), gamma=beta2s[[idx_beta2]], post_non_lin=post_non_lin, g_z=g_z, debug=debug, xz_mode=xz_mode)
    } else {
      Y <- y_from_xz(Z[,c(-1), drop=FALSE], eps_sigmaY, X=as.matrix(X_orig[,c(-1), drop=FALSE]), gamma=beta2s[[idx_beta2]], post_non_lin=post_non_lin, g_z=g_z, xz_mode=xz_mode)
    }
    is_ci <- "No_CI"
  }

  y_dir <- paste0(path_to_ukb_data, '/', is_ci, '/', n_sample[[idx_sample]], '/', seed, '/', paste0("eps_sigmaY=", eps_sigmaY))
  dir.create(y_dir, recursive = TRUE, showWarnings = FALSE)
  if (debug) {
    Y_id <- data.frame(id = X_orig$id, Y = Y$Y)
  } else {
    Y_id <- data.frame(id = X_orig$id, Y = Y)
  }
  write.csv(Y_id, file.path(y_dir, 'Y.csv'), row.names = FALSE)

  if (embedding_obs %in% c('fastsurfer', 'condVAE', 'latentDiffusion', 'freesurfer', 'medicalnet', 'cVAE_long', 'pooled_brainsynth', 'tucker_brainsynth')){
    X_obs[,c(-1)] <- scale(X_obs[,c(-1)])
    # Remove columns with NA values (zero variance columns)
    na_cols <- colSums(is.na(X_obs[, -1])) > 0
    if (any(na_cols)) {
      na_col_names <- names(X_obs[, -1])[na_cols]
      cat(sprintf("[INFO] Removing %d NA columns from X_obs: %s\n", 
                 sum(na_cols), 
                paste(na_col_names, collapse=", ")))
      X_obs <- X_obs[, c(TRUE, !na_cols)]  # Keep id column (TRUE) and non-NA columns
    }
    Y <- Y_id
  }else if(embedding_obs %in% c('scratch', 'medicalnet_ft', 'medicalnet_ft_frozen')){
    stopifnot(all.equal(X_obs$id, Y_id$id))
    X_obs$y <- Y

    # Setup paths for train results and embeddings
    embedding_dir <- paste0(y_dir, '/', embedding_obs)
    embeddings_path <- file.path(embedding_dir, '/test_embeddings.parquet')
    index_path <- file.path(embedding_dir, '/test_embeddings_index.csv')
    preds_path <- file.path(embedding_dir, '/test_predictions.parquet')
    ## Train or skip if embeddings already exist
    if (file.exists(embeddings_path)){
      cat("✅ Embedding directory already exists. Skipping training...\n")
      # Read the saved training time from previous run
      training_time_path <- file.path(embedding_dir, 'training_time.csv')
      if (file.exists(training_time_path)) {
        training_time_df <- read.csv(training_time_path)
        training_time <- training_time_df$training_time_seconds[1]
        cat(sprintf("📊 Retrieved training time from previous run: %.2f seconds\n", training_time))
      } else {
        cat("⚠️  Warning: training_time.csv not found. Setting training_time to 0.\n")
        training_time <- 0
      }
    } else {
      dir.create(embedding_dir, recursive = TRUE, showWarnings = FALSE)
      script_dir <- "inst/learn_embedding"
      # save Y and data gen config
      # Save configuration
      config <- list(
        #args = args,
        is_ci = is_ci,
        post_non_lin = post_non_lin,
        eps_sigmaX = eps_sigmaX,
        eps_sigmaY = eps_sigmaY,
        eps_sigmaZ = eps_sigmaZ,
        embedding_orig = embedding_orig,
        embedding_obs = embedding_obs,
        confounder = confounder,
        g_z = g_z,
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
        # Scale learning rate based on sample size
        base_lr_scratch <- 3e-4
        scaled_lr <- base_lr_scratch#scale_lr_by_sample_size(base_lr_scratch, n_sample[[idx_sample]], scaling_method = "power",alpha=0.05)
        
        args_vec <- c("--input_csv", normalizePath(train_csv),
                      "--id_col", "id",
                      "--output_dir", normalizePath(embedding_dir),
                      "--script_dir", normalizePath(script_dir),
                      "--model_name", "resnet18",
                      "--task", "regression",
                      "--epochs", "100",
                      "--batch_size", "16",
                      "--test_size", "0.5",
                      "--weight_decay", "0.000001",
                      "--val_frac", "0.1",
                      "--amp",
                      "--lr", sprintf("%.6e", scaled_lr),
                      "--use_tensorboard")

        start_time <- Sys.time()
        res <- run_python_safe(train_script, args_vec)
        training_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        write.csv(data.frame(training_time_seconds = training_time), file.path(embedding_dir, '/training_time.csv'), row.names = FALSE)
        cat(res, sep = "\n")  
        
      } else if (embedding_obs == 'medicalnet_ft'){
        # For 'scratch': train from scratch first, then use trained model
        cat("🏗️  Fine-tuning pretrained MedicalNet weights...\n")
          # Run training pipeline
        train_script <- "inst/learn_embedding/run_train_test_pipeline.py"
        # Scale learning rates based on sample size
        base_lr_head <- 1.2e-3
        base_lr_backbone <- 8e-5
        scaled_lr_head <- base_lr_head#scale_lr_by_sample_size(base_lr_head, n_sample[[idx_sample]], scaling_method = "log")
        scaled_lr_backbone <- base_lr_backbone#scale_lr_by_sample_size(base_lr_backbone, n_sample[[idx_sample]], scaling_method = "power", alpha = 0.05)
        
        args_vec <- c("--input_csv", normalizePath(train_csv),
                      "--id_col", "id",
                      "--output_dir", normalizePath(embedding_dir),
                      "--script_dir", normalizePath(script_dir),
                      "--model_name", "resnet18",
                      "--task", "regression",
                      "--epochs", "100",
                      "--batch_size", "16",
                      "--test_size", "0.5",
                      "--val_frac", "0.1",
                      "--amp",
                      "--lr_head", sprintf("%.6e", scaled_lr_head),
                      "--lr_backbone", sprintf("%.6e", scaled_lr_backbone),
                      "--use_tensorboard",
                      "--pretrained",
                      "--simple_head",
                      "--unfreeze_from", "layer4",
                      "--unfreeze_after_epochs", "8")

        # Use the python from the active env (auto-detected via CONDA_PREFIX)
        start_time <- Sys.time()
        res <- run_python_safe(train_script, args_vec)
        training_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        write.csv(data.frame(training_time_seconds = training_time), file.path(embedding_dir, '/training_time.csv'), row.names = FALSE)
        cat(res, sep = "\n") 
      } else if(embedding_obs == 'medicalnet_ft_frozen'){
        # For 'scratch': train from scratch first, then use trained model
        cat("🏗️  Fine-tuning pretrained MedicalNet weights with frozen backbone...\n")
          # Run training pipeline
        train_script <- "inst/learn_embedding/run_train_test_pipeline.py"
        # Scale learning rate based on sample size
        base_lr_frozen <- 1e-4
        scaled_lr <- base_lr_frozen#scale_lr_by_sample_size(base_lr_frozen, n_sample[[idx_sample]], scaling_method = "power", alpha = 0.05)
        
        args_vec <- c("--input_csv", normalizePath(train_csv),
                      "--id_col", "id",
                      "--output_dir", normalizePath(embedding_dir),
                      "--script_dir", normalizePath(script_dir),
                      "--model_name", "resnet18",
                      "--task", "regression",
                      "--epochs", "100",
                      "--batch_size", "16",
                      "--test_size", "0.5",
                      "--val_frac", "0.1",
                      "--weight_decay", "0.000001",
                      "--amp",
                      "--lr", sprintf("%.6e", scaled_lr),
                      "--use_tensorboard",
                      "--pretrained",
                      #"--simple_head",
                      "--freeze_backbone")

        # Use the python from the active env (auto-detected via CONDA_PREFIX)
        start_time <- Sys.time()
        res <- run_python_safe(train_script, args_vec)
        training_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        write.csv(data.frame(training_time_seconds = training_time), file.path(embedding_dir, '/training_time.csv'), row.names = FALSE)
        cat(res, sep = "\n") 
      }else{
        # error
        stop("Error: embedding_obs not supported")
      }
    }
    X_obs_test_trained <- as.matrix(arrow::read_parquet(embeddings_path))
    X_obs_test_trained_idx <- data.table::fread(index_path)
    X_obs <- cbind(id = X_obs_test_trained_idx$subject_id, as.data.frame(X_obs_test_trained))
    # CHECK FOR INF/NAN VALUES - COLUMN BY COLUMN:
    na_count <- sum(sapply(X_obs[,-1], function(x) sum(is.na(x))))
    inf_count <- sum(sapply(X_obs[,-1], function(x) sum(is.infinite(x) & x > 0)))
    neginf_count <- sum(sapply(X_obs[,-1], function(x) sum(is.infinite(x) & x < 0)))
    nan_count <- sum(sapply(X_obs[,-1], function(x) sum(is.nan(x))))

    cat(sprintf("[DEBUG data_gen] X_obs after loading dim=%s, NAs=%d, Inf=%d, -Inf=%d, NaN=%d\n", 
                paste(dim(X_obs), collapse="x"), na_count, inf_count, neginf_count, nan_count))

    if(inf_count > 0 || neginf_count > 0) {
      cat("[ERROR] X_obs contains Inf values after loading!\n")
      
      # Find which columns have Inf
      inf_cols <- which(sapply(X_obs[,-1], function(x) any(is.infinite(x))))
      cat(sprintf("[ERROR] Columns with Inf (indices): %s\n", paste(inf_cols, collapse=", ")))
      
      # Replace Inf with NA
      for(i in 2:ncol(X_obs)) {  # Start from 2 to skip id column
        X_obs[[i]][is.infinite(X_obs[[i]])] <- NA
      }
      cat("[FIX] Replaced Inf values with NA\n")
    }

    test_ids <- X_obs$id
    Z <- Z[match(test_ids, Z$id), ]
    Y <- Y_id[match(test_ids, Y_id$id), , drop=FALSE]
    X_orig <- X_orig[match(test_ids, X_orig$id), ]
    stopifnot(all.equal(X_obs$id, Y$id))
  }

  row.names(Z) <- 1:nrow(Z)
  row.names(X_obs) <- 1:nrow(X_obs)
  row.names(X_orig) <- 1:nrow(X_orig)
  row.names(Y) <- 1:nrow(Y)

  X_obs <- X_obs[, -c(1), drop=FALSE]
  #X_obs <- scale(X_obs)
  Z <- Z[, -c(1), drop=FALSE]
  Y <- Y[, -c(1), drop=FALSE]
  if (seed %% 100 == 0){
  cat("X_obs: ", dim(X_obs), "\n")
    cat("Y: ", dim(Y), "\n")
    cat("Z: ", dim(Z), "\n")
  }
  return(list(X_obs, Y, Z, training_time))
}


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
    X_obs <- data.table::fread(paste0(path_to_ukb_data, "/cvae_embeddings/cvae_latent_embeddings_encoder_256_mni.csv"), header=TRUE, nThread = 1)
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
    medicalnet_idx <- data.table::fread(paste0(path_to_ukb_data, "/medicalnet_embeddings/embeddings_index.csv"))
    X_obs <- cbind(id = medicalnet_idx$subject_id, as.data.frame(X_obs_emb))
  }else if(embedding_obs == 'cVAE_long'){
    X_obs <- data.table::fread(paste0(path_to_ukb_data, '/cvae_embeddings/cvae_latent_embeddings_encoder_256_penultimate.csv'), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'boc_brainsynth'){
    X_obs_emb <- arrow::read_parquet(paste0(path_to_ukb_data, '/brainsynth_embeddings/baseline_vqvae/embeddings/all_boc_embedding.parquet'))
    X_obs <- as.data.frame(X_obs_emb)
    boc_idx <- data.table::fread(paste0(path_to_ukb_data, "/brainsynth_embeddings/baseline_vqvae/embeddings/subject_ids.csv"))
    X_obs <- cbind(id = boc_idx$x, as.data.frame(X_obs_emb))
  }else if(embedding_obs == 'pooled_brainsynth'){
    X_obs_emb <- arrow::read_parquet(paste0(path_to_ukb_data, '/brainsynth_embeddings/baseline_vqvae/embeddings/all_pooled_embeddings.parquet'))
    X_obs <- as.data.frame(X_obs_emb)
    pooled_idx <- data.table::fread(paste0(path_to_ukb_data, "/brainsynth_embeddings/baseline_vqvae/embeddings/subject_ids.csv"))
    X_obs <- cbind(id = pooled_idx$x, as.data.frame(X_obs_emb))
  }else if(embedding_obs == 'tucker_brainsynth'){
    X_obs <- arrow::read_parquet(paste0(path_to_ukb_data, '/brainsynth_embeddings/baseline_vqvae/embeddings/all_tucker_embeddings.parquet'))
    colnames(X_obs)[1] <- "id"
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
  return(Z)
}

# Function to check if a column has only 2 values (binary after sd e.g.)
is_binary <- function(x) {
  if (all(is.na(x))) return(FALSE)  # Handle all-NA columns
  unique_values <- unique(x[!is.na(x)])  # Exclude NAs when checking
  length(unique_values) == 2
}

# Function to run python script safely in console
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

#' Scale learning rate based on sample size
#' 
#' @param base_lr Base learning rate (at reference sample size)
#' @param n Current sample size
#' @param n_base Reference sample size (default: 500)
#' @param scaling_method Scaling method: "sqrt" or "linear" (default: "sqrt")
#' @param alpha exponent for "power"/"capped_power"
#' @param max_factor max scaling for "capped_power"
#' @return Scaled learning rate
scale_lr_by_sample_size <- function(
  base_lr,
  n,
  n_base = 460,
  scaling_method = c("sqrt", "linear", "power", "log", "capped_power"),
  alpha = 0.3,        # exponent for "power"/"capped_power"
  max_factor = 3      # max scaling for "capped_power"
) {
  method <- match.arg(scaling_method)
  ratio  <- n / n_base

  if (ratio <= 0) stop("n and n_base must be > 0")

  scaled_lr <- switch(
    method,
    "sqrt" = base_lr * sqrt(ratio),
    "linear" = base_lr * ratio,
    "power" = base_lr * ratio^alpha,
    "log" = {
      # Gentle log scaling, factor ~ 1 + k*log(ratio)
      k <- 0.5
      base_lr * (1 + k * log(ratio))
    },
    "capped_power" = {
      raw_factor <- ratio^alpha
      factor <- min(raw_factor, max_factor)
      base_lr * factor
    }
  )

  return(scaled_lr)
}