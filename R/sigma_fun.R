# sigma_fun: updates sigma^2 using paper's closed-form expression
#'
#' Update residual variance \eqn{\sigma^2}
#'
#' Internal helper for GMERT/GMERF. Updates the residual variance parameter
#' \eqn{\sigma^2} using the closed-form expression derived from the
#' generalized log-likelihood, aggregating cluster-wise weighted residuals.
#'
#' @param N Integer. Total number of observations.
#' @param G Integer. Number of clusters.
#' @param idx List of length `G`. Row indices for each cluster.
#' @param b Numeric matrix (G x q). Current random-effects estimates.
#' @param y Numeric vector length N. Response or pseudo-response.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param D Numeric matrix (q x q). Random-effects covariance matrix.
#' @param Vi List of length `G`. Cluster-specific covariance matrices \eqn{V_i}.
#' @param s2 Numeric scalar. Current residual variance.
#' @param fhat Numeric vector length N. Fixed-effects predictions.
#' @param w Numeric vector length N. Working weights.
#'
#' @return Numeric scalar. Updated estimate of \eqn{\sigma^2}.
#'
#' @keywords internal
sigma_fun <- function(N,    # N   : total number of observations
                      G,    # G   : number of clusters
                      idx,  # idx : list of row indices per cluster
                      b,    # b   : matrix q × n vector of random effects
                      y,    # y   : response vector length N
                      Z,    # Z   : N x q random-effects design
                      D,    # D   : q x q covariance of random effects
                      Vi,   # Vi  : list of length G with V_i matrices (n_i x n_i)
                      s2,   # s2  : scalar sigma^2 (residual variance; homoskedastic)
                      fhat, # fhat: current fixed-effect fit
                      w     # length N weights
){
  term_sum <- 0
  for (g in seq_len(G)) {
    Zi <- Z[idx[[g]], , drop = FALSE]
    Wi_sq <- diag(sqrt(w[idx[[g]]]))                    # W_i^{1/2}
    e <- Wi_sq %*% y[idx[[g]]] - Wi_sq %*% fhat[idx[[g]]] - Wi_sq %*% Zi%*% b[g, ]       #  weighted residuals epsilon_i = W_i^{1/2} y - W_i^{1/2} fhat - W_i^{1/2} Z_i b_i
    term_sum <- term_sum + crossprod(e) + s2 * (length(idx[[g]]) -
                                                  s2 * sum(diag(solve(Vi[[g]]))))     # add cluster-i contribution
  }
  # Output:
  #   scalar sigma^2 new estimate
  term_sum / N                                             # divide by total N
}
