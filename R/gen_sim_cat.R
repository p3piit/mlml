#' Simulate categorical mixed-effects data for GMERT
#'
#' Generate a clustered dataset for a multinomial logit mixed-effects model with
#' group-level random intercepts and random slopes. The response has \code{K}
#' classes, with class \code{K} used as the reference category.
#'
#' @param n Integer. Total number of observations to generate.
#' @param n_groups Integer. Number of groups for group-level random effects.
#' @param n_vars Integer. Total number of predictors to simulate.
#' @param n_random Integer. Number of predictors that receive group-level random slopes.
#' @param K Integer. Number of response classes. Must be at least 3.
#' @param rho Numeric in [0, 1). Pairwise correlation used to build a
#'   compound-symmetry covariance matrix for predictors.
#' @param beta Numeric matrix or NULL. Fixed-effect coefficient matrix of dimension
#'   \code{n_vars x (K-1)}. When NULL, coefficients are sampled from N(0, 1).
#' @param mu Numeric scalar or numeric vector of length \code{n_vars}. Mean(s) used
#'   when drawing predictors.
#' @param random_slope_sd Numeric scalar or numeric vector of length \code{n_random}.
#'   Standard deviation(s) for each random slope within each contrast.
#' @param intercept_sd Numeric scalar. Standard deviation for the group intercepts
#'   within each contrast.
#' @param intercept_slope_cor Numeric scalar or numeric vector of length
#'   \code{n_random}. Correlation(s) between the intercept and each slope within
#'   each contrast.
#' @param slope_corr Numeric scalar or numeric matrix. If scalar, used as a common
#'   pairwise correlation among slopes within each contrast; if matrix, must be
#'   \code{n_random x n_random}.
#' @param D Optional full random-effects covariance matrix of dimension
#'   \code{((K-1)*(1+n_random)) x ((K-1)*(1+n_random))}. If supplied, it overrides
#'   the default block-diagonal construction.
#' @param seed Integer or NULL. Optional RNG seed for reproducibility.
#' @param random_vars Integer vector or NULL. Indices of the predictors that receive
#'   random slopes; if NULL random indices are selected.
#' @param intercept Numeric scalar or numeric vector of length \code{K-1}. Fixed
#'   intercept(s) for the non-reference classes.
#' @param return_prob Logical. If TRUE, append the true class probabilities to the
#'   output data frame.
#'
#' @details
#' For each observation \eqn{i} and non-reference class \eqn{k = 1, \dots, K-1},
#' the linear predictor is
#' \deqn{
#' \eta_{ik} = \alpha_k + x_i^\top \beta_k + z_i^\top b_{g(i)k},
#' }
#' where \eqn{g(i)} is the group of observation \eqn{i}. Class probabilities are
#' obtained through the multinomial softmax with class \eqn{K} as the reference:
#' \deqn{
#' P(Y_i = k) = \frac{\exp(\eta_{ik})}{1 + \sum_{\ell=1}^{K-1}\exp(\eta_{i\ell})},
#' \qquad
#' P(Y_i = K) = \frac{1}{1 + \sum_{\ell=1}^{K-1}\exp(\eta_{i\ell})}.
#' }
#'
#' By default, the random-effects covariance is block-diagonal across contrasts,
#' with the same within-contrast covariance repeated for each of the \code{K-1}
#' logits. To simulate cross-contrast dependence, supply a full covariance matrix
#' through \code{D}.
#'
#' @return A data.frame with columns \code{id}, \code{y}, and predictors
#'   \code{x1}, ..., \code{xP}. If \code{return_prob = TRUE}, the class
#'   probabilities \code{prob_1}, ..., \code{prob_K} are also returned.
#'

#' @importFrom MASS mvrnorm
#' @export
gen_sim_cat <- function(
  n = 1000,
  n_groups = 20,
  n_vars = 10,
  n_random = 3,
  K = 3,
  rho = 0.3,
  beta = NULL,                    # n_vars x (K-1)
  mu = 0,                         # scalar or length n_vars
  random_slope_sd = 0.5,          # scalar or length n_random
  intercept_sd = 0.8,
  intercept_slope_cor = 0.1,      # scalar or length n_random
  slope_corr = 0.1,               # scalar or n_random x n_random matrix
  D = NULL,                       # full covariance override
  seed = NULL,
  random_vars = NULL,
  intercept = 0,                  # scalar or length (K-1)
  return_prob = FALSE
) {
  if (!is.null(seed)) set.seed(seed)

  if (K < 3L) stop("K must be at least 3 for categorical simulation.")
  K1 <- K - 1L

  # interpret n_random == NULL or 0 as random-intercept only
  has_random <- !is.null(n_random) && length(n_random) == 1 && n_random > 0
  if (!has_random) {
    n_random <- 0L
  } else {
    if (n_random > n_vars) stop("n_random must be <= n_vars")
  }

  # validate random_vars when provided
  if (!is.null(random_vars)) {
    if (!all(random_vars %in% seq_len(n_vars))) stop("random_vars must be valid indices")
    if (n_random == 0L) stop("random_vars supplied but n_random is NULL/0")
    n_random <- length(random_vars)
  }

  if (rho < 0 || rho >= 1) stop("rho must be in [0, 1)")
  if (length(mu) == 1L) mu <- rep(mu, n_vars)
  if (length(mu) != n_vars) stop("mu must be a scalar or a numeric vector of length n_vars")

  # validate intercept
  if (length(intercept) == 1L) intercept <- rep(intercept, K1)
  if (length(intercept) != K1) stop("intercept must be scalar or length K-1")

  # validate random_slope_sd / intercept_slope_cor only when slopes present
  if (n_random > 0L) {
    if (length(random_slope_sd) == 1L) random_slope_sd <- rep(random_slope_sd, n_random)
    if (length(random_slope_sd) != n_random) {
      stop("random_slope_sd must be scalar or length n_random")
    }

    if (length(intercept_slope_cor) == 1L) {
      intercept_slope_cor <- rep(intercept_slope_cor, n_random)
    }
    if (length(intercept_slope_cor) != n_random) {
      stop("intercept_slope_cor must be scalar or length n_random")
    }
    if (any(intercept_slope_cor < -1 | intercept_slope_cor > 1)) {
      stop("intercept_slope_cor must be in [-1, 1]")
    }
  } else {
    random_slope_sd <- numeric(0)
    intercept_slope_cor <- numeric(0)
  }

  # group assignment (balanced as possible)
  groups <- rep(seq_len(n_groups), length.out = n)
  groups <- sample(groups, size = n, replace = FALSE)

  # predictor covariance (compound symmetry)
  Sigma_X <- matrix(rho, n_vars, n_vars)
  diag(Sigma_X) <- 1
  Xmat <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma_X)
  colnames(Xmat) <- paste0("x", seq_len(n_vars))

  # fixed-effect coefficients
  if (is.null(beta)) {
    beta <- matrix(stats::rnorm(n_vars * K1, mean = 0, sd = 1),
                   nrow = n_vars, ncol = K1)
  } else {
    beta <- as.matrix(beta)
    if (!all(dim(beta) == c(n_vars, K1))) {
      stop("beta must be an n_vars x (K-1) matrix")
    }
  }

  # choose which predictors have random slopes
  if (n_random > 0L && is.null(random_vars)) {
    random_vars <- sample(seq_len(n_vars), size = n_random)
  } else if (n_random == 0L) {
    random_vars <- integer(0)
  }

  # construct within-contrast random-effects covariance: (1 + n_random) x (1 + n_random)
  re_dim <- 1L + n_random
  Sigma_re <- matrix(0, nrow = re_dim, ncol = re_dim)
  Sigma_re[1, 1] <- intercept_sd^2

  if (n_random > 0L) {
    Sigma_re[2:re_dim, 2:re_dim] <- diag(as.matrix(random_slope_sd^2))

    # intercept-slope covariances
    Sigma_re[1, 2:re_dim] <- intercept_slope_cor * intercept_sd * random_slope_sd
    Sigma_re[2:re_dim, 1] <- Sigma_re[1, 2:re_dim]

    # slope-slope correlations
    if (length(slope_corr) == 1L) {
      if (slope_corr != 0) {
        for (i in 2:re_dim) for (j in 2:re_dim) if (i != j) {
          Sigma_re[i, j] <- slope_corr * sqrt(Sigma_re[i, i] * Sigma_re[j, j])
        }
      }
    } else {
      slope_corr <- as.matrix(slope_corr)
      if (!all(dim(slope_corr) == c(n_random, n_random))) {
        stop("slope_corr must be scalar or n_random x n_random matrix")
      }
      for (i in seq_len(n_random)) for (j in seq_len(n_random)) {
        Sigma_re[i + 1L, j + 1L] <- slope_corr[i, j] *
          random_slope_sd[i] * random_slope_sd[j]
      }
    }
  }

  # full covariance across all contrasts
  if (is.null(D)) {
    # default: independent contrasts, same Sigma_re repeated for each contrast
    D <- as.matrix(Matrix::bdiag(replicate(K1, Sigma_re, simplify = FALSE)))
  } else {
    D <- as.matrix(D)
    if (!all(dim(D) == c(K1 * re_dim, K1 * re_dim))) {
      stop("D must have dimension ((K-1)*(1+n_random)) x ((K-1)*(1+n_random))")
    }
  }

  # check positive-definiteness
  ev <- eigen(D, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) {
    stop("constructed random-effects covariance is not positive-definite; adjust SDs/correlations")
  }

  # draw group-level random effects: rows = groups, cols = stacked random effects
  re_draws <- MASS::mvrnorm(n = n_groups, mu = rep(0, K1 * re_dim), Sigma = D)

  # fixed part: N x (K-1)
  fixed_part <- Xmat %*% beta
  fixed_part <- sweep(fixed_part, 2, intercept, FUN = "+")

  # random contribution: N x (K-1)
  rand_part <- matrix(0, n, K1)

  for (i in seq_len(n)) {
    g <- groups[i]

    for (k in seq_len(K1)) {
      cols_k <- ((k - 1L) * re_dim + 1L):(k * re_dim)

      # intercept contribution for contrast k
      rand_part[i, k] <- re_draws[g, cols_k[1L]]

      # random slopes contribution for contrast k
      if (n_random > 0L) {
        rand_coefs <- re_draws[g, cols_k[-1L]]
        rand_part[i, k] <- rand_part[i, k] +
          sum(Xmat[i, random_vars] * rand_coefs)
      }
    }
  }

  # total linear predictor: N x (K-1)
  eta <- fixed_part + rand_part

  # softmax probabilities with class K as reference
  exp_eta <- exp(eta)
  denom <- 1 + rowSums(exp_eta)

  P <- matrix(0, n, K)
  P[, 1:K1] <- exp_eta / denom
  P[, K] <- 1 / denom

  # draw categorical outcome
  y_int <- apply(P, 1, function(p) sample.int(K, size = 1L, prob = p))
  y <- factor(y_int, levels = seq_len(K), labels = paste0("class", seq_len(K)))

  df <- data.frame(
    id = factor(groups),
    y = y,
    Xmat,
    stringsAsFactors = FALSE
  )

  if (return_prob) {
    colnames(P) <- paste0("prob_", seq_len(K))
    df <- cbind(df, as.data.frame(P))
  }

  df
}
