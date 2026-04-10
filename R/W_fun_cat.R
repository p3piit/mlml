#' Build cluster-specific multinomial working-weight matrices
#'
#' Construct the cluster-level working-weight matrices \eqn{W_i} for the
#' categorical GMERT pseudo-model.
#'
#' @param mu Numeric matrix of dimension \eqn{N \times K} containing the current
#'   fitted class probabilities. Each row must sum to 1.
#' @param idx List of length \eqn{G}. Each element contains the row indices of
#'   \code{mu} corresponding to one cluster.
#'
#' @return A list of length \eqn{G}. Each element is a block-diagonal matrix of
#'   dimension \eqn{((K-1)n_i) \times ((K-1)n_i)}, where \eqn{n_i} is the size
#'   of cluster \eqn{i}. The \eqn{j}-th block is the multinomial covariance
#'   matrix
#'   \deqn{\Sigma_{ij} = \mathrm{diag}(\mu_{ij}^{(1:K-1)}) -
#'         \mu_{ij}^{(1:K-1)} \mu_{ij}^{(1:K-1)\top}.}
#'
#' @details
#' For a categorical response with \eqn{K} classes and class \eqn{K} taken as
#' the reference category, the observation-level working covariance for
#' observation \eqn{(i,j)} is the \eqn{(K-1)\times(K-1)} matrix
#' \deqn{
#' [\Sigma_{ij}]_{k\ell} =
#' \begin{cases}
#' \mu^{(k)}_{ij}(1-\mu^{(k)}_{ij}), & k=\ell,\\
#' -\mu^{(k)}_{ij}\mu^{(\ell)}_{ij}, & k\neq\ell.
#' \end{cases}
#' }
#' The cluster-specific weight matrix \eqn{W_i} is obtained by stacking these
#' blocks along the diagonal.
#'
#' This is the direct multivariate analogue of the binary working weights
#' \eqn{w_{ij} = \mu_{ij}(1-\mu_{ij})}.
#'
#' @keywords internal
W_fun_cat <- function(mu,   # mu  : N x K matrix of class probabilities
                      idx,  # idx : list length G with row indices for each cluster
                      eps = 1e-6  # ridge regularization parameter
) {
  if (!is.matrix(mu)) {
    mu <- as.matrix(mu)
  }

  N <- nrow(mu)
  K <- ncol(mu)
  K1 <- K - 1L

  if (K < 2L) {
    stop("mu must have at least 2 columns (K >= 2).")
  }

  if (any(!is.finite(mu))) {
    stop("mu contains non-finite values.")
  }

  if (any(mu < 0 | mu > 1)) {
    stop("All entries of mu must lie in [0, 1].")
  }

  row_sums <- rowSums(mu)
  if (any(abs(row_sums - 1) > 1e-8)) {
    stop("Each row of mu must sum to 1.")
  }

  W_list <- vector("list", length(idx))

  for (g in seq_along(idx)) {
    idg <- idx[[g]]
    ni <- length(idg)

    blocks <- vector("list", ni)

    for (j in seq_len(ni)) {
      p <- mu[idg[j], 1:K1, drop = TRUE]   # keep only non-reference classes

      # Sigma_ij = diag(p) - p p^T
      Sigma_ij <- diag(p, nrow = K1, ncol = K1) - tcrossprod(p)

      # --- SAFEGUARD: ridge regularization ---
      Sigma_ij <- Sigma_ij + eps * diag(K1)

      blocks[[j]] <- Sigma_ij
    }

    if (ni == 1L) {
      W_list[[g]] <- blocks[[1L]]
    } else {
      W_list[[g]] <- as.matrix(Matrix::bdiag(blocks))
    }
  }

  W_list
}