#' Simulate response Y from Z and optionally X
#'
#' @param X (non-observed) feature representation used to simulate response
#' @param Z observed confounders
#' @param eps_sigmaY integer of (homogeneous) standard deviation of noise in Y
#' @param post_non_lin integer specifying nonlinear function g() applied to input consisting of X, Z and epsY
#'
#' @return returns function to construct Y
#' @export
y_from_xz <- function(Z, eps_sigmaY,X=NULL, post_non_lin = 1) {
  g <- post_non_lin_g(post_non_lin)
  epsY <- stats::rnorm(nrow(Z), 0, eps_sigmaY)
  if(is.null(X)){
    Y <- g(scale(rowMeans(as.matrix(Z)))+epsY)
  } else {
    epsX <- stats::rnorm(nrow(Z), 0, 1)
    X <- scale(X)
    Z <- scale(Z)
    Y <- g(scale(rowMeans(as.matrix(cbind(X,Z)))+epsY))
  }
}

post_non_lin_g <- function(post_non_lin){
  if (post_non_lin %in% 1:3) {
    g <- function(s) {
      y <- s^post_non_lin
      return(y)
    }
  } else if (post_non_lin %in% 4:5) {
    g <- function(s) {
      s <- scale(s)
      y <- exp(-s^2/2) * sin(ifelse(post_non_lin == 4, 3, 24) * s)
      return(y)
    }
  }
  return(g)
}
