# sigma_fun_cat: updates residual variance sigma2 using precomputed A_i^{-1}
#'
#' Update residual variance using precomputed A_i^{-1}
#'
#' @param N Integer. Total number of observations.
#' @param G Integer. Number of clusters.
#' @param idx List of length G. Each element contains the row indices for a given
#'   cluster.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param W List of length G. Each element is the cluster-specific working weight
#'   matrix W_i of dimension ((K-1)n_i) x ((K-1)n_i).
#' @param y_t Numeric matrix (N x K1). Stacked pseudo-response stored by contrast.
#' @param fhat Numeric matrix (N x K1). Current fitted fixed-part predictions.
#' @param b Numeric matrix (G x (K1*q)). Current random effects.
#' @param Ainv List of length G. Each element is a ((K-1)q x (K-1)q) matrix
#'   representing A_i^{-1} for the corresponding cluster.
#' @param K1 Integer. Number of log-odds contrasts, i.e. K - 1.
#'
#' @return A numeric scalar giving the updated residual variance estimate.
#'
#' @keywords internal
sigma_fun_cat <- function(N,
                                G,
                                idx,
                                Z,
                                W,
                                y_t,     # N x K1
                                fhat,    # N x K1
                                b,       # G x (K1*q)
                                Ainv,
                                K1) {

  q <- ncol(Z)
  rss_total <- 0

  for (g in seq_len(G)) {
    idg <- idx[[g]]
    Zi  <- Z[idg, , drop = FALSE]                  # n_i x q
    Zti <- kronecker(diag(K1), Zi)                # ((K-1)n_i) x ((K-1)q)
    Wi  <- W[[g]]

    resid_mat <- y_t[idg, , drop = FALSE] - fhat[idg, , drop = FALSE]
    resid_vec <- c(resid_mat)                     # stack by column

    ri <- resid_vec - as.vector(Zti %*% b[g, ])   # ((K-1)n_i) x 1

    rss_total <- rss_total +
      as.numeric(crossprod(ri, Wi %*% ri)) +
      sum(diag(Ainv[[g]] %*% crossprod(Zti, Wi %*% Zti)))
  }

   max(rss_total / (N * K1), 1e-6)  # enforce a minimum variance to avoid numerical issues
}