#' Simulate response Y from Z and optionally X
#'
#' @param X (non-observed) feature representation used to simulate response
#' @param Z observed confounders
#' @param eps_sigmaY integer of (homogeneous) standard deviation of noise in Y
#' @param post_non_lin integer specifying nonlinear function g() applied to input consisting of X, Z and epsY
#'
#' @return returns function to construct Y
#' @export
y_from_xz <- function(Z, eps_sigmaY,X=NULL,beta2s=NULL, idx_beta2=NULL, post_non_lin = 1, g_z = "linear") {
  g <- post_non_lin_g(post_non_lin)
  epsY <- stats::rnorm(nrow(Z), 0, eps_sigmaY)

  # transform Z and aggregate over columns with random effects
  transformed_z <- transform_z(as.data.frame(Z), g_z)
  beta_Z_ <- stats::rnorm(ncol(transformed_z), 0, 1)
  active_Z <- stats::rbinom(ncol(transformed_z), 1, 1)
  #if(sum(active_Z)==0){
  #  active_Z[1] <- 1
  #}
  beta_Z <- abs(beta_Z_*active_Z)
  beta_Z <- beta_Z/sum(beta_Z)
  predictor_z <- as.matrix(transformed_z)%*%beta_Z

  if(is.null(X)){
    Y <- g(predictor_z+epsY)
  } else {
    beta_X_ <- stats::rnorm(ncol(X), 0, 1)
    active_X <- stats::rbinom(ncol(X), 1, 0.5)
    beta_X <- abs(beta_X_*active_X)
    beta_X <- beta_X/sum(beta_X)
    #beta_X <- rep(1/ncol(X), ncol(X))
    predictor_x <- X%*%beta_X*beta2s[[idx_beta2]]

    Y <- g(predictor_z+predictor_x+epsY)
  }
  return(Y)
}

post_non_lin_g <- function(post_non_lin){
  if (post_non_lin %in% 1:3) {
    g <- function(s) {
      y <- s^post_non_lin
      return(y)
    }
  } else if (post_non_lin == 4) {
    g <- function(s) {
      y <- tanh(s)
      return(y)
    }
  } else if (post_non_lin == 5) {
    g <- function(s) {
      y <- exp(-abs(s))
      return(y)
    }
  }
  return(g)
}

transform_z <- function(Z, g_z){
  if (g_z %in% "linear") {
    return(Z)
  } else if (g_z %in% "squared") {
    #square numeric columns
    for (col in names(Z)) {
      if (!is_binary(Z[[col]])) {
        # Add squared column if not binary
        Z[[paste0(col, "_squared")]] <- Z[[col]]^2
      }
    }
    return(Z)
  } else if (g_z %in% "realistic") {
    #square numeric columns
    for (col in names(Z)) {
      if (!is_binary(Z[[col]])) {
        # Add squared column if not binary
        Z[[paste0(col, "_squared")]] <- Z[[col]]^2
      }
    }
    # Add interaction terms for all pairs of confounders
    if ('sex' %in% names(Z)){
      Z <- fastDummies::dummy_cols(Z, select_columns='sex', remove_first_dummy = TRUE)
      if (!is_binary(Z[[col]])) {
        Z[[paste0(col, "_sex")]] <- Z[['sex_1']] * Z[[col]]
      }
      Z[['sex_1']] <- NULL
    }
    # add higher order terms for date_diff
    if("date_diff" %in% names(Z)){
      Z[['date_diff_cubed']] <- Z[['date_diff^3']]
      Z[['date_diff_fourth']] <- Z[['date_diff^4']]
    }
    return(Z)
  } else if (g_z %in% "breakpoint3"){
    for (col in names(Z)) {
      if (!is_binary(Z[[col]])) {
        # Add squared column if not binary
        Z[[paste0(col, "exp_sin")]] <- exp(-Z[[col]]^2)*sin(3*Z[[col]])
      }
    }
    return(Z)
  }
}

# Function to check if a column has only 2 values (binary after sd e.g.)
is_binary <- function(x) {
  unique_values <- unique(x)
  length(unique_values) == 2
}
