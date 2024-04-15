#' Title
#'
#' @param s either Z+epsilon or X+Z+epsilon
#'
#' @return returns function to construct Y
#' @export
y_from_xz <- function(s, post_non_lin = 1) {
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
