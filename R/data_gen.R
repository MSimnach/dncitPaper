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
#' @param response response variable
#'
#' @return list of X_obs, Y, Z
#' @export
data_gen <- function(seed, idx_sample=NULL, n_sample=NULL, idx_beta2=NULL, beta2s=NULL, n=NULL, post_non_lin=4, eps_sigmaX=0, eps_sigmaY=1, eps_sigmaZ=0, embedding_orig='fastsurfer',
                     embedding_obs='fastsurfer', confounder='AS', response='simulated'){
  set.seed(seed)

  X_orig <- load_X_orig(embedding_orig)
  X_obs <- load_X_obs(embedding_obs, embedding_orig, X_orig, eps_sigmaX)
  Z <- load_Z(confounder)

  # take only common IDs of X, Z
  colnames(X_orig)[1] <- 'id'
  colnames(X_obs)[1] <- 'id'
  colnames(Z)[1] <- 'id'
  X_orig <- as.data.frame(X_orig)
  X_obs <- as.data.frame(X_obs)
  Z <- as.data.frame(Z)
  X_merged <- dplyr::inner_join(X_orig, X_obs, by='id')
  subset_X_obs <- X_obs[X_obs$id %in% X_merged$id, ]
  subset_X_orig <-  X_orig[X_orig$id %in% X_merged$id, ]
  subset_Z <- Z[Z$id %in% X_merged$id, ]

  #sample rows
  if(is.null(beta2s)){
    idx <- sample(1:nrow(X_merged), n_sample[[idx_sample]])
  }else if(is.null(n_sample)){
    idx <- sample(1:nrow(X_merged), n)
  }
  X_obs <- subset_X_obs[idx,-c(1)]
  X_orig <- subset_X_orig[idx,-c(1)]
  Z <- subset_Z[idx,-c(1)]

  epsZ <- matrix(stats::rnorm((nrow(Z)*ncol(Z)), 0, eps_sigmaZ),nrow=nrow(Z),ncol=ncol(Z))
  Z <- Z+epsZ

  #Standardize
  Z <- scale(Z)
  X_orig <- scale(X_orig)
  X_obs <- scale(X_obs)

  if(response=='simulated'){
    if(is.null(beta2s)){
      Y <- y_from_xz(Z, eps_sigmaY, post_non_lin)
    }else if(is.null(n_sample)){
      Y <- y_from_xz(Z, eps_sigmaY, X=X_orig,post_non_lin=post_non_lin)
    }

  }

  return(list(X_obs,Y,Z))
}


load_X_orig <- function(embedding_orig){
  if(embedding_orig == 'fastsurfer'){
    X <- data.table::fread(paste("./Data/ukb_fastsurfer.csv"), header=TRUE, nThread = 1)
  }else if(embedding_orig == 'condVAE'){
    X <- data.table::fread(paste("./Data/ukb_condVAE.csv"), header=TRUE, nThread = 1)
  }else if(embedding_orig == 'latentDiffusion'){
    X <- data.table::fread(paste("./Data/ukb_latentDiffusion.csv"), header=TRUE, nThread = 1)
  }
  return(X)
}

load_X_obs <- function(embedding_obs, embedding_orig, X_orig, eps_sigmaX){
  if(embedding_obs==embedding_orig){
    X_obs <- X_orig
  }else if(embedding_obs == 'fastsurfer'){
    X_obs <- data.table::fread(paste("./Data/ukb_fastsurfer.csv"), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'condVAE'){
    X_obs <- data.table::fread(paste("./Data/ukb_condVAE.csv"), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'latentDiffusion'){
    X_obs <- data.table::fread(paste("./Data/ukb_latentDiffusion.csv"), header=TRUE, nThread = 1)
  }else if(grepl('noisy',embedding_obs, fixed=TRUE)){
    epsX <- stats::rnorm(nrow(X_orig)*(ncol(X_orig)-1), 0,eps_sigmaX)
    X_obs <- cbind(X_orig[,1], X_orig[,2:ncol(X_orig)]+epsX)
  }
}

load_Z <- function(confounder){
  if(confounder=='AS'){
    Z <- data.table::fread(paste("./Data/ukb_Z_age_sex.csv"), header=TRUE)
  }else if(confounder == 'genes10'){
    Z <- data.table::fread(paste("./Data/ukb_Z_genes10.csv"), header=TRUE)
  }
}
