#' Update cluster-level random effects
#'
#' Internal update step for GMERT/GMERF. Computes the cluster-specific random
#' effects \eqn{b_g} given current covariance components and working quantities.
#'
#' @param G Integer. Number of clusters.
#' @param Z Numeric matrix (N x q). Random-effects design matrix (stacked by cluster).
#' @param Vi List of length `G`. Each element is an \eqn{n_g \times n_g} covariance
#'   matrix for cluster `g`.
#' @param D Numeric matrix (q x q). Random-effects covariance matrix.
#' @param idx List of length `G`. Row indices for each cluster in `Z` / `y_t` / `fhat` / `w`.
#' @param y_t Numeric vector length N. Working (pseudo-)response.
#' @param fhat Numeric vector length N. Fixed-part prediction.
#' @param w Numeric vector length N. Working weights (typically \eqn{\mu (1-\mu)}).
#'
#' @return A numeric matrix (G x q) of updated random effects.
#'
#' @keywords internal
b_fun <- function(G, Z, Vi, D, idx, y_t, fhat, w) {
  b <- matrix(0, G, ncol(Z))
  for (g in seq_len(G)) {
    Zi   <- Z[idx[[g]], , drop = FALSE]
    Wi_sq <- diag(sqrt(w[idx[[g]]]))

    Zi_w <- Wi_sq %*% Zi
    rhs  <- Wi_sq %*% y_t[idx[[g]]] - Wi_sq %*% fhat[idx[[g]]]

    b[g, ] <- D %*% t(Zi_w) %*% solve(Vi[[g]]) %*% rhs
  }
  b
}
