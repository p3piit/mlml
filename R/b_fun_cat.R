##############################################################
# b_fun_cat: updates random effects b_i using precomputed A_i^{-1}
# for categorical GMERT
#'
#' Update cluster-specific random effects using precomputed A_i^{-1}
#'
#' Computes the updated random effects b_i for each cluster i using
#' the precomputed inverse matrices A_i^{-1} for the categorical GMERT.
#'
#' @param G Integer. Number of clusters.
#' @param Z Numeric matrix (N x q). Random-effects design matrix for one contrast.
#' @param W List of length G. Each element contains the cluster-specific working
#'   weight matrix W_i of dimension ((K-1)n_i) x ((K-1)n_i).
#' @param idx List of length G. Each element contains the row indices of Z
#'   corresponding to a given cluster.
#' @param y_t Numeric vector of length N*(K-1). Stacked pseudo-response.
#' @param fhat Numeric vector of length N*(K-1). Stacked fitted fixed-part predictions.
#' @param Ainv List of length G. Each element is a ((K-1)q x (K-1)q) matrix
#'   representing A_i^{-1} for the corresponding cluster.
#' @param sigma2 Numeric scalar. Current residual variance.
#' @param K1 Integer. Number of log-odds contrasts, i.e. K - 1.
#'
#' @return A numeric matrix of dimension (G x ((K-1)q)) containing the updated
#'   random effects for each cluster.
#'
#' @details
#' For each cluster i, this function computes
#'   b_i = A_i^{-1} (1/sigma2) Z_tilde_i^T W_i (y_tilde_i - f_i),
#' where Z_tilde_i = I_(K-1) \otimes Z_i.
#'
#' @keywords internal
b_fun_cat <- function(G,
                            Z,        # N x q random-effects design matrix
                            W,       # list of length G with cluster-specific weight matrices W_i
                            idx,     # list of length G with row indices for each cluster
                            y_t,      # N x K1 matrix
                            fhat,     # N x K1 matrix
                            Ainv,    # list of length G with A_i^{-1} matrices
                            sigma2,  # residual variance
                            K1) {    # number of contrasts = K - 1
  q <- ncol(Z)
  b <- matrix(0, G, K1 * q)

  for (g in seq_len(G)) {
    idg <- idx[[g]]

    Zi  <- Z[idg, , drop = FALSE]                # n_i x q
    Zti <- kronecker(diag(K1), Zi)              # ((K-1)n_i) x ((K-1)q)
    Wi  <- W[[g]]

    resid_mat <- y_t[idg, , drop = FALSE] - fhat[idg, , drop = FALSE]  # n_i x K1
    resid <- as.vector(as.matrix(resid_mat))    # stacks by column: class 1, ..., class K-1

    rhs <- (1 / sigma2) * crossprod(Zti, Wi %*% resid)
    b[g, ] <- as.vector(Ainv[[g]] %*% rhs)
  }

  b
}
