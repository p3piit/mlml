# Vi_fun: builds list of V_i = Z_i D Z_i^T + sigma^2 I_{n_i}
#'
#' Build cluster-specific covariance matrices V_i
#'
#' Internal helper for GMERT/GMERF. For each cluster \eqn{i}, constructs
#' \eqn{V_i = (W_i^{1/2} Z_i) D (W_i^{1/2} Z_i)^T + \sigma^2 I_{n_i}},
#' where \eqn{W_i} is the diagonal working-weight matrix.
#'
#' @param Z Numeric matrix (N x q). Random-effects design matrix (stacked by cluster).
#' @param D Numeric matrix (q x q). Random-effects covariance matrix.
#' @param s2 Numeric scalar. Residual variance on the working scale.
#' @param G Integer. Number of clusters.
#' @param idx List of length `G`. Row indices for each cluster.
#' @param w Numeric vector length N. Working weights.
#'
#' @return A list of length `G`. Each element is an \eqn{n_i \times n_i}
#'   covariance matrix \eqn{V_i}.
#'
#' @keywords internal
Vi_fun <- function(Z,   # Z   : N x q random-effects design
                   D,   # D   : q x q covariance of random effects
                   s2,  # s2  : scalar sigma^2 (residual variance; homoskedastic)
                   G,   # G   : number of clusters
                   idx, # idx : list of row indices per cluster
                   w    # length N weights
){
  Vi <- list()
  for (g in seq_len(G)) {
    Zi <- Z[idx[[g]], , drop = FALSE]                                    # Z_i
    Wi_sq <- diag(sqrt(w[idx[[g]]]))                                     # W_i^{1/2}
    Zi_w <- Wi_sq %*% Zi                                                 # W^{1/2} Z_i
    Vi[[g]] <- Zi_w %*% D %*% t(Zi_w) + as.numeric(s2) * diag(nrow(Zi_w))  # V_i
  }
  # Output:
  # list of length G with n_i x n_i covariance matrices V_i
  return(Vi)
}
