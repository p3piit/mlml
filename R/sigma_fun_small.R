# sigma_fun_Ainv: updates residual variance sigma2 using precomputed A_i^{-1}
#'
#' Update residual variance using precomputed \eqn{A_i^{-1}}
#'
#' Updates the residual variance parameter \eqn{\sigma^2} using cluster-wise
#' quantities and the precomputed matrices \eqn{A_i^{-1}} (typically obtained
#' from \code{\link{Ajnv_fun}}). This is the “small” (optimized) counterpart of
#' \code{\link{sigma_fun}}, avoiding direct formation/inversion of \eqn{V_i}.
#'
#' @param N Integer. Total number of observations.
#' @param G Integer. Number of clusters.
#' @param idx List of length G. Each element contains the row indices for a given
#'   cluster.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param w Numeric vector of length N. Observation-level weights.
#' @param y_t Numeric vector of length N. Pseudo-response (working response).
#' @param fhat Numeric vector of length N. Current fitted fixed-part predictions.
#' @param b Numeric matrix (G x q). Current random effects (one row per cluster).
#' @param Ainv List of length G. Each element is a (q x q) matrix representing
#'   \eqn{A_i^{-1}} for the corresponding cluster.
#' @param sigma2 Numeric scalar. Current residual variance value used in the update.
#'
#' @return A numeric scalar giving the updated residual variance estimate.
#'
#' @details
#' The update aggregates a weighted residual sum of squares across clusters and
#' includes a trace adjustment term involving \eqn{\mathrm{tr}(V_i^{-1})},
#' computed efficiently via \eqn{A_i^{-1}}.
#'
#' @keywords internal
sigma_fun_small <- function(N,        # N      : total number of rows
                            G,        # G      : number of clusters
                            idx,      # idx    : list length G with row indices for each cluster
                            Z,        # Z      : N x q random-effects design
                            w,        # W      : list of weigths
                            y_t,      # y_t    : pseudo-response vector
                            fhat,     # fhat   : fitted fixed-part prediction
                            b,        # b      : matrix G x q with current random effects
                            Ainv,     # Ainv   : list length G; each Ainv[[g]] = A_i^{-1}
                            sigma2   # sigma2 : residual variance (current)
){
  rss_total <- 0

  for (g in seq_len(G)) {
    Zi  <- Z[idx[[g]], , drop = FALSE]             # n_i x q
    wi  <- diag(w[idx[[g]]])                             # weights matrix
    ni  <- length(wi)

    # residuals eps_i
    eps_i <- diag(sqrt(w[idx[[g]]])) %*% (y_t[idx[[g]]] - fhat[idx[[g]]] - as.vector(Zi %*% b[g, ]))

    # tr(V_i^{-1})
    trVi_inv <- (ni / sigma2) - as.numeric((1 / sigma2^2)) * sum(diag(Ainv[[g]] %*% t(Zi) %*% wi %*% Zi))

    # accumulate rss
    rss_total <- rss_total + t(eps_i) %*% eps_i + sigma2 * (ni - sigma2 * trVi_inv)
  }

  rss_total / N
}
