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
                     embedding_obs='fastsurfer', confounder='AS', g_z='linear', debug=FALSE, xz_mode='independent'){
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
  Z[,c(-1)] <- Z[,c(-1)]+epsZ

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

  X_orig[,c(-1)] <- scale(X_orig[,c(-1)])
  #X_obs <- scale(X_obs)

  if(is.null(beta2s)){
    if (debug) {
      Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z, debug=debug)
    } else {
      Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z)
    }
    is_ci <- "CI"
  }else{
    if (debug) {
      Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, X=as.matrix(X_orig[,c(-1)]), gamma=beta2s[[idx_beta2]], post_non_lin=post_non_lin, g_z=g_z, debug=debug, xz_mode=xz_mode)
    } else {
      Y <- y_from_xz(Z[,c(-1)], eps_sigmaY, X=as.matrix(X_orig[,c(-1)]), gamma=beta2s[[idx_beta2]], post_non_lin=post_non_lin, g_z=g_z, xz_mode=xz_mode)
    }
    is_ci <- "No_CI"
  }

  y_dir <- paste0(path_to_ukb_data, '/', is_ci, '/', n_sample[[idx_sample]], '/', seed, '/', paste0("eps_sigmaY=", eps_sigmaY))
  dir.create(y_dir, recursive = TRUE)
  if (debug) {
    Y_id <- data.frame(id = X_orig$id, Y = Y$Y)
  } else {
    Y_id <- data.frame(id = X_orig$id, Y = Y)
  }
  write.csv(Y_id, file.path(y_dir, 'Y.csv'), row.names = FALSE)

  if (embedding_obs %in% c('fastsurfer', 'condVAE', 'latentDiffusion', 'freesurfer', 'medicalnet')){
    X_obs[,c(-1)] <- scale(X_obs[,c(-1)])
  }else if(embedding_obs %in% c('scratch', 'medicalnet_ft', 'medicalnet_ft_frozen')){
    stopifnot(all.equal(X_obs$id, Y_id$id))
    X_obs$y <- Y

    # Setup paths for train results and embeddings
    embedding_dir <- paste0(y_dir, '/', embedding_obs)
    embeddings_path <- file.path(embedding_dir, '/test_embeddings.parquet')
    index_path <- file.path(embedding_dir, '/test_embeddings_index.csv')
    preds_path <- file.path(embedding_dir, '/test_predictions.parquet')
    ## Train or skip if embeddings already exist
    if (file.exists(embedding_dir)){
      cat("✅ Embedding directory already exists. Skipping training...\n")
    } else {
      dir.create(embedding_dir, recursive = TRUE)
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

  X_obs <- X_obs[, -c(1)]
  X_obs <- scale(X_obs)
  Z <- Z[, -c(1)]
  Y <- Y[, -c(1)]
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
    medicalnet_idx <- fread(paste0(out_dir, "/embeddings_index.csv"))
    X_obs <- cbind(id = medicalnet_idx$subject_id, as.data.frame(X_obs_emb))
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
  unique_values <- unique(x)
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