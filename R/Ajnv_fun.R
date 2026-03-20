#' Precompute cluster-specific inverse matrices for random-effects updates
#'
#' Computes the inverse matrices \eqn{A_j^{-1}} for each cluster \eqn{j}, where
#' \eqn{A_j = D^{-1} + (1 / \sigma^2) Z_j^T W_j Z_j}. These matrices are used to
#' speed up the closed-form updates of random effects and variance components
#' in the small-sample GMERT / GMERF algorithms.
#'
#' @param G Integer. Number of clusters.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param w Numeric vector of length N. Observation-level weights.
#' @param D Numeric matrix (q x q). Current random-effects covariance matrix.
#' @param sigma2 Numeric scalar. Current residual variance.
#' @param idx List of length G. Each element contains the row indices of \code{Z}
#'   corresponding to a given cluster.
#'
#' @return A list of length G, where each element is a (q x q) matrix containing
#'   \eqn{A_j^{-1}} for the corresponding cluster.
#'
#' @details
#' The matrices returned by this function are used in alternative update formulas
#' for the random effects \eqn{b_j}, the residual variance \eqn{\sigma^2}, and the
#' random-effects covariance matrix \eqn{D}, reducing computational cost by
#' avoiding repeated matrix inversions.
#'
#' @keywords internal
Ajnv_fun <- function(G,        # G      : number of clusters
                     Z,        # Z      : N x q random-effects design
                     w,        # W      : list of weigths
                     D,        # D      : q x q covariance of random effects (current)
                     sigma2,   # sigma2 : residual variance (current)
                     idx      # idx    : list length G with row indices for each cluster
){
  Ajnv_list <- vector("list", G)

  Djnv <- solve(D)                                 # D^{-1} (q x q) once

  for (g in seq_len(G)) {
    Zj  <- Z[idx[[g]], , drop = FALSE]             # n_j x q
    if (length(w[idx[[g]]]) ==1) {
      Wj <- as.matrix(w[idx[[g]]])
    } else {
      Wj  <- diag(w[idx[[g]]])
    }
    ZtWZ <- t(Zj) %*% Wj %*% Zj                     # Z_j^T W_j Z_j (q x q)

    A <- Djnv + as.numeric((1 / sigma2)) * ZtWZ                # A_j
    Ajnv <- solve(A)                              # A_j^{-1}

    Ajnv_list[[g]] <- Ajnv
  }

  Ajnv_list
}
