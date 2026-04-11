#' Precompute cluster-specific inverse matrices for random-effects updates
#'
#' Computes the inverse matrices A_j^{-1} for each cluster j in the
#' categorical GMERT, where
#'   A_j = D^{-1} + (1 / sigma2) * Z_tilde_j' W_j Z_tilde_j.
#'
#' @param G Integer. Number of clusters.
#' @param Z Numeric matrix (N x q). Random-effects design matrix for one contrast.
#' @param W List of length G. Each element is the cluster-specific working weight
#'   matrix W_j of dimension ((K-1)*n_j) x ((K-1)*n_j).
#' @param D Numeric matrix ((K-1)*q x (K-1)*q). Current random-effects covariance matrix.
#' @param sigma2 Numeric scalar. Current residual variance.
#' @param idx List of length G. Each element contains the row indices of Z
#'   corresponding to a given cluster.
#' @param K1 Integer. Number of log-odds contrasts, i.e. K - 1.
#'
#' @return A list of length G, where each element is a ((K-1)*q x (K-1)*q)
#'   matrix containing A_j^{-1} for the corresponding cluster.
#'
#' @keywords internal
Ajnv_fun_cat <- function(G,        # G      : number of clusters
                         Z,        # Z      : N x q random-effects design
                         W,        # W      : list of cluster-specific weight matrices
                         D,        # D      : ((K-1)q x (K-1)q) covariance of random effects
                         sigma2,   # sigma2 : residual variance
                         idx,      # idx    : list length G with row indices for each cluster
                         K1        # K1     : number of contrasts = K - 1
) {
  Ajnv_list <- vector("list", G)

  Djnv <- solve(D)  # D^{-1} once

  for (g in seq_len(G)) {
    Zj <- Z[idx[[g]], , drop = FALSE]           # n_j x q
    Zt_j <- kronecker(diag(K1), Zj)             # Z_tilde_j = I_(K-1) \otimes Z_j
    Wj <- W[[g]]                                # ((K-1)n_j) x ((K-1)n_j)

    A <- Djnv + (1 / sigma2) * crossprod(Zt_j, Wj %*% Zt_j)
    Ajnv_list[[g]] <- solve(A)
  }

  Ajnv_list
}