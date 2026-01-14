# gll_fun: computes generalized log-likelihood (up to additive constant)
#'
#' Compute generalized log-likelihood (up to an additive constant)
#'
#' Internal helper for GMERT/GMERF. Computes the generalized log-likelihood
#' criterion used for monitoring and convergence checks. The returned value is
#' defined up to an additive constant (i.e., it is suitable for comparing
#' iterations but not as an absolute likelihood value).
#'
#' @param idx List of row indices per cluster.
#' @param b Numeric matrix (G x q). Random-effects estimates.
#' @param y Numeric vector length N. Response or pseudo-response.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param D Numeric matrix (q x q). Random-effects covariance matrix.
#' @param s2 Numeric scalar. Residual variance.
#' @param fhat Numeric vector length N. Fixed-effects predictions.
#' @param w Numeric vector length N. Working weights.
#'
#' @return Numeric scalar. Generalized log-likelihood value (up to an additive constant).
#'
#' @keywords internal
gll_fun <- function(idx, # idx : list of row indices per cluster
                    b,   # b   : matrix q × n vector of random effects
                    y,   # y   : response vector length N
                    Z,   # Z   : N x q random-effects design
                    D,   # D   : q x q covariance of random effects
                    s2,  # s2  : scalar sigma^2 (residual variance; homoskedastic)
                    fhat,# fhat: current fixed-effect fit
                    w    # length N weights
){
  logdetD <- as.numeric(determinant(D, logarithm = TRUE)$modulus)     # log|D|
  term_b <- sum(rowSums((b %*% solve(D)) * b))                         # sum_i b_i^T D^{-1} b_i
  term_r <- 0
  term_logRi <- 0
  for (i in seq_along(idx)) {
    Zi <- Z[idx[[i]], , drop = FALSE]                      # Z_i
    Wi_sq <- diag(sqrt(w[idx[[i]]]))                       # W_i^{1/2}
    e <- Wi_sq %*% y[idx[[i]]] - Wi_sq %*% fhat[idx[[i]]] - Wi_sq %*% Zi%*% b[i, ]       #  weighted residuals epsilon_i = W_i^{1/2} y - W_i^{1/2} fhat - W_i^{1/2} Z_i b_i
    term_r <- term_r + sum(e^2) / s2                                  # (1/s2) * ||ri||^2
    term_logRi <- term_logRi + length(idx[[i]]) * log(s2)              # n_i * log(s2)
  }

  # Output:
  #   scalar GLL value (up to constant)
  # GLL = sum_i [ (y - f - Zb)^T R_i^{-1} (y - f - Zb) + log|R_i| ] +
  # + sum_i b_i^T D^{-1} b_i + G*log|D|
  term_r + term_logRi + term_b + length(idx) * logdetD                 # assemble GLL
}
