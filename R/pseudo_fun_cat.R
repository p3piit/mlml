#' Build multinomial pseudo-response for categorical GMERT
#'
#' Construct the observation-level pseudo-response matrix used in the
#' PQL / IRLS update for the categorical GMERT.
#'
#' @param Y Numeric matrix of dimension \eqn{N \times K} containing the one-hot
#'   encoded response. Each row must contain exactly one 1.
#' @param mu Numeric matrix of dimension \eqn{N \times K} containing the current
#'   fitted class probabilities. Each row must sum to 1.
#' @param eta Numeric matrix of dimension \eqn{N \times (K-1)} containing the
#'   current linear predictors for the non-reference classes.
#' @param eps Numeric scalar. Ridge safeguard added to the diagonal of each
#'   multinomial covariance block before inversion. Default is 1e-5.
#' @param N Integer. Number of observations.
#' @param K1 Integer. Number of contrasts = K - 1.
#'
#' @return A numeric matrix of dimension \eqn{N \times (K-1)} containing the
#'   pseudo-response for the non-reference classes.
#'
#' @details
#' For each observation \eqn{i}, let \eqn{y_i^*}, \eqn{\mu_i^*}, and
#' \eqn{\eta_i} denote the \eqn{(K-1)}-vectors obtained by removing the
#' reference class. The pseudo-response is computed as
#' \deqn{
#' \tilde y_i = \eta_i + \Sigma_i^{-1}(y_i^* - \mu_i^*)
#' }
#' where
#' \deqn{
#' \Sigma_i = \mathrm{diag}(\mu_i^*) - \mu_i^* {\mu_i^*}^\top.
#' }
#'
#' A ridge term \eqn{\varepsilon I} is added to \eqn{\Sigma_i} before inversion
#' for numerical stability.
#'
#' This is the multivariate analogue of the binary pseudo-response
#' \eqn{\eta + (y-\mu)/[\mu(1-\mu)]}.
#'
#' @keywords internal
pseudo_fun_cat <- function(Y,      # Y   : N x K one-hot response matrix
                           mu,     # mu  : N x K fitted class probabilities
                           eta,    # eta : N x (K-1) current linear predictor
                           N,      # N   : number of observations
                           K1,     # K1 = K-1 number of contrasts
                           eps = 1e-5
) {
  y_t <- matrix(0, N, K1)

  for (i in seq_len(N)) {
    yi  <- Y[i, 1:K1, drop = TRUE]     # observed non-reference classes
    mui <- mu[i, 1:K1, drop = TRUE]    # fitted probs, non-reference classes
    etai <- eta[i, , drop = TRUE]      # current linear predictor

    Sigma_i <- diag(mui, nrow = K1, ncol = K1) - tcrossprod(mui)

    # safeguard for inversion
    Sigma_i <- Sigma_i + eps * diag(K1)

    y_t[i, ] <- etai + solve(Sigma_i, yi - mui)
  }

  y_t
}