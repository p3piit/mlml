# D_fun: updates D using the paper's expression (average of b_i b_i^T plus shrinkage term)
#'
#' Update random-effects covariance matrix D
#'
#' Internal helper for GMERT/GMERF. Updates the random-effects covariance
#' matrix \eqn{D} using the closed-form expression based on the empirical
#' second moments of the random effects and a shrinkage correction term
#' derived from the conditional covariance.
#'
#' @param G Integer. Number of clusters.
#' @param idx List of length `G`. Row indices for each cluster.
#' @param b Numeric matrix (G x q). Current random-effects estimates.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param D Numeric matrix (q x q). Current random-effects covariance matrix.
#' @param Vi List of length `G`. Cluster-specific covariance matrices \eqn{V_i}.
#' @param w Numeric vector length N. Working weights.
#'
#' @return Numeric matrix (q x q). Updated estimate of the random-effects
#'   covariance matrix \eqn{D}.
#'
#' @keywords internal
D_fun <- function(G,   # G   : number of clusters
                  idx, # idx : list of row indices per cluster
                  b,   # b   : matrix q × n vector of random effects
                  Z,   # Z   : N x q random-effects design
                  D,   # D   : q x q covariance of random effects
                  Vi,  # Vi  : list of length G with V_i matrices (n_i x n_i)
                  w    # length N weights
){
  q <- ncol(Z)
  S <- matrix(0, q, q)                                     # accumulator
  for (g in seq_len(G)) {
    Zi <- Z[idx[[g]], , drop = FALSE]                      # Z_i
    Wi_sq <- diag(sqrt(w[idx[[g]]]))                       # W_i^{1/2}
    Zi_w <- Wi_sq %*% Zi                                   # W^{1/2} Z_i
    S <- S + tcrossprod(b[g, ]) +                          # b_i b_i^T
      (D - D %*% t(Zi_w) %*% solve(Vi[[g]]) %*% Zi_w %*% D)    # shrinkage term per cluster
  }
  # Output:
  #   q x q updated D
  S/G                                                      # average over clusters
}
