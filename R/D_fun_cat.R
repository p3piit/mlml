# D_fun_cat: updates random-effects covariance D using precomputed A_i^{-1}
#'
#' Update random-effects covariance matrix using precomputed A_i^{-1}
#'
#' @param G Integer. Number of clusters.
#' @param b Numeric matrix (G x ((K-1)q)). Current random effects,
#'   one row per cluster.
#' @param Ainv List of length G. Each element is a (((K-1)q) x ((K-1)q)) matrix
#'   representing A_i^{-1} for the corresponding cluster.
#'
#' @return A numeric (((K-1)q) x ((K-1)q)) matrix giving the updated estimate
#'   of the random-effects covariance matrix D.
#'
#' @keywords internal
D_fun_cat <- function(G,
                            b,
                            Ainv) {
  q_star <- ncol(b)
  D_new <- matrix(0, q_star, q_star)

  for (g in seq_len(G)) {
    D_new <- D_new + tcrossprod(b[g, ]) + Ainv[[g]]
  }

  D_new / G
}