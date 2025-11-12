#' Simulate response Y from confounders Z and optionally features X
#'
#' Generates a response variable Y based on confounders Z and optionally features X.
#' The function supports multiple modes for controlling the relationship between X and Z effects,
#' allows for nonlinear transformations, and provides flexible control over signal strength and mixing.
#'
#' @param Z Matrix or data.frame of observed confounders (n x q)
#' @param eps_sigmaY Numeric. Standard deviation of homogeneous noise in Y
#' @param X Matrix of feature representations (n x p). If NULL, Y depends only on Z
#' @param post_non_lin Integer (1-5). Specifies nonlinear function g() applied to linear predictor:
#'   \itemize{
#'     \item 1: identity (s)
#'     \item 2: quadratic (s^2)
#'     \item 3: cubic (s^3)
#'     \item 4: tanh(s)
#'     \item 5: exp(-|s|)
#'   }
#' @param g_z Character. Transformation applied to Z before computing predictor:
#'   "linear", "squared", "realistic", or "breakpoint3"
#' @param xz_mode Character. Controls relationship between X and Z effects:
#'   \itemize{
#'     \item "independent": beta_X and beta_Z are independent
#'     \item "orthogonalize": X is orthogonalized w.r.t. Z before drawing beta_X
#'     \item "correlate": enforces correlation between beta_X and beta_Z via rho
#'   }
#' @param gamma Numeric in [0,1]. Mixing weight between Z and X effects:
#'   Y ~ g((1-gamma)*pZ + gamma*pX + eps)
#' @param tau Numeric. Overall signal strength scaling factor
#' @param sparsity_Z Numeric in [0,1]. Probability of non-zero coefficients in beta_Z
#' @param sparsity_X Numeric in [0,1]. Probability of non-zero coefficients in beta_X
#' @param shrink_X Logical. Whether to apply shrinkage to correlation matrix of X
#' @param lambda_X Numeric. Shrinkage parameter for correlation matrix (if shrink_X=TRUE)
#' @param rho Numeric in [-1,1]. Target correlation between beta_X and beta_Z (used when xz_mode="correlate")
#' @param debug Logical. If TRUE, returns list with Y and diagnostic info; if FALSE, returns only Y vector
#'
#' @return If debug=FALSE: numeric vector Y of length n.
#'   If debug=TRUE: list with components:
#'   \itemize{
#'     \item Y: numeric vector of responses
#'     \item info: list containing pZ, pX, beta_Z, beta_X, correlation, and parameters
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simulate Y from Z only (conditional independence)
#' Z <- matrix(rnorm(100 * 5), 100, 5)
#' Y <- y_from_xz(Z, eps_sigmaY = 1, post_non_lin = 4)
#'
#' # Simulate Y from both X and Z with gamma mixing
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' Y <- y_from_xz(Z, eps_sigmaY = 1, X = X, gamma = 0.3,
#'                xz_mode = "independent", post_non_lin = 4)
#'
#' # Get diagnostic information
#' result <- y_from_xz(Z, eps_sigmaY = 1, X = X, gamma = 0.5, debug = TRUE)
#' cor(result$info$pX, result$info$pZ)
#' }
y_from_xz <- function(Z, eps_sigmaY, X = NULL,
                      post_non_lin = 1, g_z = "linear",
                      xz_mode = "Sigma=I_p",
                      gamma = 0.5, tau = 1.0,
                      sparsity_Z = 1.0, sparsity_X = 0.2,
                      shrink_X = TRUE, lambda_X = 0.05,
                      rho = 0.0, debug = FALSE) {
                        
  g <- post_non_lin_g(post_non_lin)
  n <- nrow(Z)
  epsY <- stats::rnorm(n, 0, eps_sigmaY)

  # Z branch
  Z_star <- transform_z(as.data.frame(Z), g_z)
  Z_star <- scale_numeric_df(Z_star)
  Z_mat  <- as.matrix(Z_star)

  beta_Z <- stats::rnorm(ncol(Z_mat), 0, 1)
  if (sparsity_Z < 1) {
    beta_Z <- beta_Z * stats::rbinom(length(beta_Z), 1, sparsity_Z)
  }
  pZ_raw <- drop(Z_mat %*% beta_Z)
  pZ     <- as.numeric(scale(pZ_raw))

  # If no X provided, done
  if (is.null(X)) {
    s <- tau * pZ
    Y <- g(s + epsY)
    if (debug) {
      return(list(Y = Y, info = list(pZ = pZ, beta_Z = beta_Z, mode = "Z-only")))
    } else {
      return(Y)
    }
  }

  # if X is provided
  X_mat <- as.matrix(X)
  if (xz_mode == "orthogonalize") {
    # Remove linear overlap of X with Z
    X_use <- residualize_X_on_Z(X_mat, Z_mat)
    beta_X <- sample_beta_corr(X_use, sparsity = sparsity_X, shrink = shrink_X, lambda = lambda_X)
    pX_raw <- drop(X_use %*% beta_X)

  } else if (xz_mode == "correlate") {
    # Build a mapping W: Z -> X by OLS (min ||X - ZW||_F^2)
    # W = (Z'Z)^(-1) Z' X  (use QR solve)
    W <- qr.solve(Z_mat, X_mat)                # q x p
    mapped <- drop(t(W) %*% beta_Z)            # p-vector in X-space aligned with beta_Z
    # Normalize components and combine with rho
    mapped <- mapped / sqrt(sum(mapped^2) + 1e-12)
    indep  <- sample_beta_corr(X_mat, sparsity = sparsity_X, shrink = shrink_X, lambda = lambda_X)
    indep  <- indep / sqrt(sum(indep^2) + 1e-12)
    rho <- max(-1, min(1, rho))
    beta_X <- rho * mapped + sqrt(1 - rho^2) * indep
    # rescale to original magnitude ~ N(0, Cor(X))
    beta_X <- beta_X * sd(sample_beta_corr(X_mat, sparsity = 1.0, shrink = shrink_X, lambda = lambda_X))
    pX_raw <- drop(X_mat %*% beta_X)

  } else if (xz_mode == "independent") { # "independent"
    beta_X <- sample_beta_corr(X_mat, sparsity = sparsity_X, shrink = shrink_X, lambda = lambda_X)
    pX_raw <- drop(X_mat %*% beta_X)
  } else { # "Sigma=I_p", i.e. beta_X with Sigma=I_p
    beta_X_ <- stats::rnorm(ncol(X_mat), 0, 1)
    active_X <- stats::rbinom(ncol(X_mat), 1, sparsity_X)
    beta_X <- beta_X_*active_X
    pX_raw <- drop(X_mat %*% beta_X)
  }

  pX <- as.numeric(scale(pX_raw))

  # Mix balanced predictors
  gamma <- max(0, min(1, gamma))
  s <- tau * ((1 - gamma) * pZ + gamma * pX)

  Y <- g(s + epsY)

  info <- list(
    pZ = pZ, pX = pX,
    beta_Z = beta_Z, beta_X = beta_X,
    mu_Z = mean(pZ), sd_Z = sd(pZ),
    mu_X = mean(pX), sd_X = sd(pX),
    corr_pX_pZ = suppressWarnings(stats::cor(pX, pZ)),
    mode = xz_mode, gamma = gamma, tau = tau, rho = rho
  )
  if (debug) {
    return(list(Y = Y, info = info))
  } else {
    return(Y)
  }
}

# ==============================
# Helpers
# ==============================

# Tiny PD-safe Cholesky with optional shrinkage
safe_chol <- function(S, shrink = TRUE, lambda = 0.05) {
  if (shrink) {
    S <- (1 - lambda) * S + lambda * diag(ncol(S))
  }
  # try Cholesky; if it fails, go nearPD (Matrix) and retry
  tryCatch(
    chol(S),
    error = function(e) {
      S_pd <- as.matrix(Matrix::nearPD(S, corr = TRUE)$mat)
      chol(S_pd)
    }
  )
}

# Sample beta ~ N(0, Cor(X)) with sparsity mask
sample_beta_corr <- function(X, sparsity = 1.0, shrink = TRUE, lambda = 0.05) {
  stopifnot(is.matrix(X))
  p <- ncol(X)
  S <- stats::cor(X, use = "complete.obs")
  L <- t(safe_chol(S, shrink = shrink, lambda = lambda))  # lower-tri by transposing upper-tri
  z <- stats::rnorm(p)
  beta <- drop(L %*% z)
  if (sparsity < 1) {
    mask <- stats::rbinom(p, 1, sparsity)
    beta <- beta * mask
  }
  beta
}

# Orthogonalize X wrt Z (QR projection)
residualize_X_on_Z <- function(X, Z) {
  # Z: n x q (full rank), X: n x p
  qrZ <- qr(Z)
  Q  <- qr.Q(qrZ)
  X - Q %*% (t(Q) %*% X)
}

# Scale numeric columns of a data.frame (leave binary as-is)
scale_numeric_df <- function(Df) {
  Df2 <- Df
  for (nm in names(Df2)) {
    x <- Df2[[nm]]
    if (!is_binary(x)) Df2[[nm]] <- as.numeric(scale(x))
  }
  Df2
}

# Treat a vector as binary if it has exactly two distinct values
is_binary <- function(x) {
  length(unique(x)) == 2
}

# Nonlinear post-transform
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

# Confounder transforms (yours, with small fixes)
transform_z <- function(Z, g_z){
  Z <- as.data.frame(Z)
  if (g_z == "linear") {
    return(Z)
  } else if (g_z == "squared") {
    for (col in names(Z)) {
      if (!is_binary(Z[[col]])) Z[[paste0(col, "_squared")]] <- Z[[col]]^2
    }
    return(Z)
  } else if (g_z == "realistic") {
    for (col in names(Z)) {
      if (!is_binary(Z[[col]])) Z[[paste0(col, "_squared")]] <- Z[[col]]^2
    }
    if ('sex' %in% names(Z)){
      if (!requireNamespace("fastDummies", quietly = TRUE)) {
        stop("Please install 'fastDummies' for realistic sex dummies.")
      }
      Z <- fastDummies::dummy_cols(Z, select_columns='sex', remove_first_dummy = TRUE)
      sexcol <- 'sex_1'
      for (col in names(Z)) {
        if (col != sexcol && !is_binary(Z[[col]]) && sexcol %in% names(Z)) {
          Z[[paste0(col, "_sex")]] <- Z[[sexcol]] * Z[[col]]
        }
      }
      if (sexcol %in% names(Z)) Z[[sexcol]] <- NULL
    }
    if ("date_diff" %in% names(Z)){
      Z[['date_diff_cubed']]  <- Z[['date_diff']]^3
      Z[['date_diff_fourth']] <- Z[['date_diff']]^4
    }
    return(Z)
  } else if (g_z == "breakpoint3"){
    for (col in names(Z)) {
      if (!is_binary(Z[[col]])) {
        Z[[paste0(col, "exp_sin")]] <- exp(-Z[[col]]^2) * sin(3*Z[[col]])
      }
    }
    return(Z)
  }
  Z
}
