# gll_fun_cat: computes generalized log-likelihood (up to additive constant)
#'
#' Compute generalized log-likelihood (up to an additive constant)
#'
#' Internal helper for categorical GMERT/GMERF. Computes the Gaussian
#' pseudo-model criterion used for monitoring and convergence checks.
#'
#' @param idx List of row indices per cluster.
#' @param b Numeric matrix (G x (K1*q)). Random-effects estimates.
#' @param y Numeric matrix (N x K1). Stacked pseudo-response by contrast.
#' @param Z Numeric matrix (N x q). Random-effects design matrix.
#' @param D Numeric matrix ((K1*q) x (K1*q)). Random-effects covariance matrix.
#' @param s2 Numeric scalar. Residual variance.
#' @param fhat Numeric matrix (N x K1). Fixed-effects predictions.
#' @param W List of length G. Cluster-specific working weight matrices.
#' @param K1 Integer. Number of contrasts = K - 1.
#'
#' @return Numeric scalar. Generalized log-likelihood value
#'   (up to an additive constant).
#'
#' @keywords internal
gll_fun_cat <- function(idx,
                        b,
                        y,      # N x K1 pseudo-response
                        Z,
                        D,
                        s2,
                        fhat,   # N x K1 fitted fixed part
                        W,
                        K1) {

  q <- ncol(Z)
  G <- length(idx)

  D_inv <- solve(D)
  logdetD <- as.numeric(determinant(D, logarithm = TRUE)$modulus)

  term_b <- sum(rowSums((b %*% D_inv) * b))
  term_r <- 0
  term_logRi <- 0

  for (g in seq_along(idx)) {
    idg <- idx[[g]]
    ni  <- length(idg)

    Zi  <- Z[idg, , drop = FALSE]
    Zti <- kronecker(diag(K1), Zi)             # ((K-1)n_i) x ((K-1)q)
    Wi  <- W[[g]]

    y_vec    <- c(y[idg, , drop = FALSE])      # stack by column
    fhat_vec <- c(fhat[idg, , drop = FALSE])

    r_i <- y_vec - fhat_vec - as.vector(Zti %*% b[g, ])

    term_r <- term_r + as.numeric(crossprod(r_i, Wi %*% r_i)) / s2
    term_logRi <- term_logRi + (K1 * ni) * log(s2)
    # optional fuller version: term_logRi <- term_logRi +
    #   (K1 * ni) * log(s2) - as.numeric(determinant(Wi, logarithm = TRUE)$modulus)
  }

  term_r + term_logRi + term_b + G * logdetD
}