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
                     embedding_obs='fastsurfer', confounder='AS', g_z='linear'){
  #path_to_ukb_data <- ""
  set.seed(seed)

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
  X_obs <- subset_X_obs[idx,-c(1)]
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

  row.names(Z) <- 1:nrow(Z)
  row.names(X_obs) <- 1:nrow(X_obs)
  row.names(X_orig) <- 1:nrow(X_orig)
  row.names(Y) <- 1:nrow(Y)

  return(list(X_obs,Y,Z))
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
  }
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


