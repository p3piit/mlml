##############################################################
# b_fun_Ainv: updates random effects b_i using precomputed A_i^{-1}
#'
#' Update cluster-specific random effects using precomputed \eqn{A_i^{-1}}
#'
#' Computes the updated random effects \eqn{b_i} for each cluster \eqn{i} using
#' the precomputed inverse matrices \eqn{A_i^{-1}} (typically returned by
#' \code{\link{Ajnv_fun}}). This version avoids forming and inverting \eqn{V_i}
#' directly, and is used in the “small” GMERT / GMERF optimization routines.
#'
#' @param G Integer. Number of clusters.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param w Numeric vector of length N. Observation-level weights.
#' @param idx List of length G. Each element contains the row indices of \code{Z}
#'   (and associated vectors) corresponding to a given cluster.
#' @param y_t Numeric vector of length N. Pseudo-response (working response).
#' @param fhat Numeric vector of length N. Current fitted fixed-part predictions.
#' @param Ainv List of length G. Each element is a (q x q) matrix representing
#'   \eqn{A_i^{-1}} for the corresponding cluster.
#' @param sigma2 Numeric scalar. Current residual variance.
#'
#' @return A numeric matrix of dimension (G x q) containing the updated random
#'   effects for each cluster.
#'
#' @details
#' For each cluster \eqn{i}, this function computes
#' \deqn{b_i = A_i^{-1} \cdot (1/\sigma^2) Z_i^T W_i (y^*_i - f_i)}
#' where \eqn{y^*_i} is the pseudo-response restricted to cluster \eqn{i},
#' \eqn{f_i} is the fixed-part prediction, and \eqn{W_i} is a diagonal matrix of
#' working weights.
#'
#' @keywords internal
b_fun_small <- function(G,        # G      : number of clusters
                        Z,        # Z      : N x q random-effects design
                        w,        # W      : list of weigths
                        idx,      # idx    : list length G with row indices for each cluster
                        y_t,      # y_t    : pseudo-response vector (length N)
                        fhat,     # fhat   : fitted fixed-part prediction (length N)
                        Ainv,     # Ainv   : list length G; each Ainv[[g]] = A_i^{-1} (q x q)
                        sigma2    # sigma2 : residual variance (current)
){
  q <- ncol(Z)
  b <- matrix(0, G, q)

  for (g in seq_len(G)) {
    Zi   <- Z[idx[[g]], , drop = FALSE]           # n_i x q
    Wi   <- diag(w[idx[[g]]])                     # weights vector (length n_i)
    resid<- y_t[idx[[g]]] - fhat[idx[[g]]]        # (y_t - fhat) on cluster i

    rhs  <- as.numeric((1 / sigma2)) * t(Zi) %*% (Wi %*% resid)  # right-hand side
    b[g, ] <- as.vector(Ainv[[g]] %*% rhs)           # b_i = A_i^{-1} * rhs
  }

  b
}
