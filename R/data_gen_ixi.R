## Data generation for IXI dataset with neural network embeddings
library(dplyr)
library(data.table)
#' IXI Data Generation with Neural Network Embeddings
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
data_gen_ixi <- function(seed, idx_sample=NULL, n_sample=NULL, idx_beta2=NULL, beta2s=NULL, 
                         post_non_lin=4, eps_sigmaX=0, eps_sigmaY=1, eps_sigmaZ=0, 
                         embedding_orig='medicalnet', embedding_obs='scratch', confounder='ukb_z4', 
                         g_z='linear'){
  
  path_to_ixi_data <- Sys.getenv("IXI_PATH", unset = NA)
  
  if (is.na(path_to_ixi_data) || path_to_ixi_data == "") {
    stop("❌ Environment variable IXI_PATH is not set. Please define it in your .Renviron.")
  }
  
  set.seed(seed)
  
  # Load original data (CSV with image paths and y values)
  X_orig <- load_X_orig_ixi(path_to_ixi_data, embedding_orig)
  
  # Load observed data (neural network embeddings)
  X_obs <- load_X_obs_ixi(path_to_ixi_data, embedding_obs, embedding_orig, X_orig, eps_sigmaX)
  
  # Load confounders
  Z <- load_Z_ixi(path_to_ixi_data, confounder)
  
  # take only common IDs of X, Z
  colnames(X_orig)[1] <- 'id'
  colnames(X_obs)[1] <- 'id'
  colnames(Z)[1] <- 'id'
  
  # Convert IXI IDs to numeric for consistency
  numeric_ids <- as.integer(gsub("IXI", "", X_orig$id))
  X_orig$id <- numeric_ids
  numeric_ids <- as.integer(gsub("IXI", "", X_obs$id))
  X_obs$id <- numeric_ids
  
  X_orig <- as.data.frame(X_orig)
  X_obs <- as.data.frame(X_obs)
  Z <- as.data.frame(Z)
  
  # Remove duplicated IDs in Z
  Z <- Z[!duplicated(Z$id), ]
  
  common_ids <- Reduce(intersect, list(X_orig$id, X_obs$id, Z$id))
  
  subset_X_obs <- X_obs[X_obs$id %in% common_ids, ]
  subset_X_orig <- X_orig[X_orig$id %in% common_ids, ]
  subset_Z <- Z[Z$id %in% common_ids, ]
  
  # Check if ids are equal
  stopifnot(all.equal(subset_X_orig[,1], subset_Z[,1]))
  
  # Sample rows
  idx <- sample(1:nrow(subset_Z), n_sample[[idx_sample]])
  X_obs <- subset_X_obs[idx, -c(1)]
  X_orig <- subset_X_orig[idx, -c(1)]
  Z <- subset_Z[idx, -c(1), drop=FALSE]
  
  # Add noise to confounders
  epsZ <- matrix(stats::rnorm((nrow(Z)*ncol(Z)), 0, eps_sigmaZ), nrow=nrow(Z), ncol=ncol(Z))
  Z <- Z + epsZ
  
  # Remove zero columns and multicollinearity in one-hot-encoding of sites (due to subsampling)
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
  # Remove one site column
  if(site_sum == nrow(Z)){
    for(site in site_columns){
      if(is_binary(Z[[site]])){
        Z <- Z[, !names(Z) %in% site]
        break
      }
    }
  }
  
  # Standardize
  # Scale only continuous confounders
  Z <- Z %>% dplyr::mutate(dplyr::across(dplyr::where(function(x) !is_binary(x)), scale))
  
  # For X_orig, we use the y values from the CSV (already continuous)
  X_orig <- scale(X_orig)
  
  # For X_obs, we use the neural network embeddings (already continuous)
  X_obs <- scale(X_obs)
  
  # Generate Y from confounders and potentially X_orig
  if(is.null(beta2s)){
    Y <- y_from_xz(Z, eps_sigmaY, post_non_lin=post_non_lin, g_z=g_z)
  } else {
    Y <- y_from_xz(Z, eps_sigmaY, X=X_orig, beta2s=beta2s, idx_beta2=idx_beta2, post_non_lin=post_non_lin, g_z=g_z)
  }
  
  row.names(Z) <- 1:nrow(Z)
  row.names(X_obs) <- 1:nrow(X_obs)
  row.names(X_orig) <- 1:nrow(X_orig)
  row.names(Y) <- 1:nrow(Y)
  
  return(list(X_obs, Y, Z))
}


load_X_orig_ixi <- function(path_to_ixi_data, embedding_orig){
  if(embedding_orig == 'medicalnet'){
    # Load the CSV with image paths and y values
    X <- data.table::fread(paste0(path_to_ixi_data, '/t1/csv/all_labeled.csv'), header=TRUE, nThread = 1)
    # For X_orig, we use the y values as the original features
    X_orig <- data.frame(id = X$id, y = X$y)
    return(X_orig)
  }
  # Add other embedding types if needed
  return(X)
}


load_X_obs_ixi <- function(path_to_ixi_data, embedding_obs, embedding_orig, X_orig, eps_sigmaX){
  if(embedding_obs %in% c('medicalnet', 'scratch', 'ft_head_only', 'ft_full')){
    
    # Check if embeddings already exist
    embedding_dir <- paste0(path_to_ixi_data, '/t1/embeddings/', embedding_obs)
    embeddings_path <- file.path(embedding_dir, 'embeddings.npy')
    index_path <- file.path(embedding_dir, 'embeddings_index.csv')
    
    if(file.exists(embeddings_path) && file.exists(index_path)){
      cat("📂 Loading existing embeddings from:", embeddings_path, "\n")
      
      # Load embeddings using reticulate
      library(reticulate)
      np <- import("numpy")
      
      # Load embeddings and index
      embeddings <- np$load(embeddings_path)
      index_df <- read.csv(index_path, stringsAsFactors = FALSE)
      
      # Convert to data frame
      embedding_df <- as.data.frame(embeddings)
      colnames(embedding_df) <- paste0("embed_", 1:ncol(embedding_df))
      
      # Add IDs
      X_obs <- data.frame(id = index_df$subject_id, embedding_df)
      
      cat("✅ Loaded embeddings with shape:", dim(embeddings), "\n")
      
    } else {
      cat("🔄 Extracting new embeddings using neural network...\n")
      
      # Create embedding directory
      dir.create(embedding_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Use reticulate to run the Python embedding extraction script
      library(reticulate)
      
      # Get the path to the extraction script
      script_path <- file.path(dirname(dirname(path_to_ixi_data)), 
                              "Coding/dncitPaper/inst/learn_embedding/extract_ixi_embeddings.py")
      
      # Input CSV path
      input_csv <- paste0(path_to_ixi_data, '/t1/csv/all_labeled.csv')
      
      # Determine model name based on embedding type
      model_name <- switch(embedding_obs,
                          'medicalnet' = 'resnet18',
                          'scratch' = 'resnet18', 
                          'ft_head_only' = 'resnet18',
                          'ft_full' = 'resnet18',
                          'resnet18')  # default
      
      # Run the Python script
      cmd <- sprintf("python %s --input_csv %s --output_dir %s --model_name %s --batch_size 4",
                     script_path, input_csv, embedding_dir, model_name)
      
      cat("Running command:", cmd, "\n")
      result <- system(cmd, intern = TRUE)
      
      if(length(result) > 0){
        cat("Python script output:\n")
        cat(paste(result, collapse = "\n"), "\n")
      }
      
      # Check if extraction was successful
      if(file.exists(embeddings_path) && file.exists(index_path)){
        cat("✅ Embedding extraction completed successfully\n")
        
        # Load the newly created embeddings
        np <- import("numpy")
        embeddings <- np$load(embeddings_path)
        index_df <- read.csv(index_path, stringsAsFactors = FALSE)
        
        # Convert to data frame
        embedding_df <- as.data.frame(embeddings)
        colnames(embedding_df) <- paste0("embed_", 1:ncol(embedding_df))
        
        # Add IDs
        X_obs <- data.frame(id = index_df$subject_id, embedding_df)
        
        cat("✅ Loaded new embeddings with shape:", dim(embeddings), "\n")
        
      } else {
        stop("❌ Failed to extract embeddings. Check Python environment and script.")
      }
    }
    
    return(X_obs)
    
  } else if(grepl('noisy', embedding_obs, fixed=TRUE)){
    # Add noise to original features
    epsX <- stats::rnorm(nrow(X_orig)*(ncol(X_orig)-1), 0, eps_sigmaX)
    X_obs <- cbind(X_orig[,1], scale(X_orig[,2:ncol(X_orig)])+epsX)
    return(X_obs)
  }
  
  # Default: return original
  return(X_orig)
}


load_Z_ixi <- function(path_to_ixi_data, confounder){
  if(confounder == 'ukb_z4'){
    # Load IXI confounder data
    Z <- data.table::fread(paste0(path_to_ixi_data, "/ixi_confounder.csv"), header=TRUE, nThread = 1)
    return(Z)
  }
  # Add other confounder types if needed
  stop("Confounder type not supported for IXI data: ", confounder)
}


# Function to check if a column has only 2 values (binary after sd e.g.)
is_binary <- function(x) {
  unique_values <- unique(x)
  length(unique_values) == 2
}
