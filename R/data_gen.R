## Data generation
#' Title
#'
#' @param seed random seed for current data generation
#' @param idx_sample index of current sample size
#' @param n_sample list of sample sizes
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
data_gen <- function(seed, idx_sample, n_sample, post_non_lin=4, eps_sigmaX=0, eps_sigmaY=1, eps_sigmaZ=0, embedding_orig='fastsurfer',
                     embedding_obs='fastsurfer', confounder='AS', response='simulated'){
  set.seed(seed)
  if(post_non_lin==1){
    g = function(s){
      return(s)
    }
  }else if(post_non_lin==2){
    g = function(s){
      return(s^2)
    }
  }else if(post_non_lin==3){
    g = function(s){
      return(s^3)
    }
  }else if(post_non_lin==4){
    g = function(s){
      s = scale(s)
      temp = exp(-s^2/2)*sin(3*s)
      return(temp)
    }
  }else if(post_non_lin==5){
    g = function(s){
      s = scale(s)
      temp = exp(-s^2/2)*sin(24*s)
      return(temp)
    }
  }
  if(embedding_orig == 'fastsurfer'){
    X_orig <- data.table::fread(paste("./Data/ukb_fastsurfer.csv"), header=TRUE, nThread = 1)
  }else if(embedding_orig == 'condVAE'){
    X_orig <- data.table::fread(paste("./Data/ukb_condVAE.csv"), header=TRUE, nThread = 1)
  }else if(embedding_orig == 'latentDiffusion'){
    X_orig <- data.table::fread(paste("./Data/ukb_latentDiffusion.csv"), header=TRUE, nThread = 1)
  }

  if(embedding_obs==embedding_orig){
    X_obs <- X_orig
  }else if(embedding_obs == 'fastsurfer'){
    X_obs <- data.table::fread(paste("./Data/ukb_fastsurfer.csv"), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'condVAE'){
    X_obs <- data.table::fread(paste("./Data/ukb_condVAE.csv"), header=TRUE, nThread = 1)
  }else if(embedding_obs == 'latentDiffusion'){
    X_obs <- data.table::fread(paste("./Data/ukb_latentDiffusion.csv"), header=TRUE, nThread = 1)
  }else if(grepl('noisy',embedding_obs, fixed=TRUE)){
    epsX <- rnorm(nrow(X_orig)*(ncol(X_orig)-1), 0,eps_sigmaX)
    X_obs <- cbind(X_orig[,1], X_orig[,2:ncol(X_orig)]+epsX)
  }

  if(confounder=='AS'){
    Z <- data.table::fread(paste("./Data/ukb_Z_age_sex.csv"), header=TRUE)
  }else if(confounder == 'genes10'){
    Z <- data.table::fread(paste("./Data/ukb_Z_genes10.csv"), header=TRUE)
  }

  # take only common IDs of X, Z
  colnames(X_orig)[1] <- 'id'
  colnames(X_obs)[1] <- 'id'
  colnames(Z)[1] <- 'id'
  X_orig <- as.data.frame(X_orig)
  X_obs <- as.data.frame(X_obs)
  Z <- as.data.frame(Z)
  X_merged <- inner_join(X_orig, X_obs, by='id')
  subset_X_obs <- X_obs[X_obs$id %in% X_merged$id, ]
  subset_X_orig <-  X_orig[X_orig$id %in% X_merged$id, ]
  subset_Z <- Z[Z$id %in% X_merged$id, ]

  #sample rows
  idx <- sample(1:nrow(X_merged), n_sample[[idx_sample]])
  X_obs <- subset_X_obs[idx,-c(1)]
  X_orig <- subset_X_orig[idx,-c(1)]
  Z <- subset_Z[idx,-c(1)]

  epsZ <- matrix(rnorm((nrow(Z)*ncol(Z)), 0, eps_sigmaZ),nrow=nrow(Z),ncol=ncol(Z))
  Z <- Z+epsZ

  Z <- scale(Z)
  X_orig <- scale(X_orig)
  X_obs <- scale(X_obs)

  if(response=='simulated'){
    epsY <- rnorm(nrow(Z), 0, eps_sigmaY)
    Y <- g(scale(rowMeans(as.matrix(Z)))+epsY)
  }

  return(list(X_obs,Y,Z))
}
