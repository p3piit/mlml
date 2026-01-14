# D_fun_small: updates random-effects covariance D using precomputed A_i^{-1}
#'
#' Update random-effects covariance matrix using precomputed \eqn{A_i^{-1}}
#'
#' Updates the random-effects covariance matrix \eqn{D} by combining the outer
#' products of the current random effects \eqn{b_i} with the corresponding
#' precomputed inverse matrices \eqn{A_i^{-1}}. This function is the optimized
#' (“small”) counterpart of \code{D_fun}, avoiding direct dependence on \eqn{V_i}.
#'
#' @param G Integer. Number of clusters.
#' @param b Numeric matrix (G x q). Current random effects, one row per cluster.
#' @param Ainv List of length G. Each element is a (q x q) matrix representing
#'   \eqn{A_i^{-1}} for the corresponding cluster.
#'
#' @return A numeric (q x q) matrix giving the updated estimate of the
#'   random-effects covariance matrix \eqn{D}.
#'
#' @details
#' The update is computed as
#' \deqn{D = \frac{1}{G} \sum_{i=1}^G \left( b_i b_i^T + A_i^{-1} \right),}
#' which combines the empirical covariance of the random effects with a
#' shrinkage term derived from the inverse information matrices.
#'
#' @keywords internal
D_fun_small <- function(G,        # G      : number of clusters
                        b,        # b      : matrix G x q with current random effects
                        Ainv      # Ainv   : list length G; each Ainv[[g]] = A_i^{-1}
){
  q <- ncol(b)
  D_new <- matrix(0, q, q)

  for (g in seq_len(G)) {
    D_new <- D_new + tcrossprod(b[g, ]) + Ainv[[g]]   # b_i b_i^T + A_i^{-1}
  }

  D_new / G
}
