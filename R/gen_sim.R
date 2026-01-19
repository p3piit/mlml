#' Simulate mixed-effects data with correlated random effects
#'
#' Generate a dataset for a linear mixed model with group-level random
#' intercepts and random slopes. Random effects for the intercept and the
#' slopes are drawn from a multivariate normal whose covariance matrix is
#' constructed from supplied standard deviations and correlations.
#'
#' @param n Integer. Total number of observations to generate.
#' @param n_groups Integer. Number of groups for group-level random effects.
#' @param n_vars Integer. Total number of predictors to simulate.
#' @param n_random Integer. Number of predictors that receive group-level random slopes.
#' @param rho Numeric in [0, 1). Pairwise correlation used to build a compound‑symmetry
#'   covariance matrix for predictors (controls multicollinearity).
#' @param beta Numeric vector or NULL. Fixed-effect coefficients of length `n_vars`.
#'   When NULL coefficients are sampled from N(0, 1).
#' @param mu Numeric scalar or numeric vector length `n_vars`. Mean(s) used when
#'   drawing predictors; a scalar is recycled to length `n_vars`.
#' @param random_slope_sd Numeric scalar or numeric vector length `n_random`.
#'   Standard deviation(s) for each random slope.
#' @param intercept_sd Numeric. Standard deviation for the group intercepts.
#' @param intercept_slope_cor Numeric scalar or numeric vector length `n_random`.
#'   Correlation(s) between the intercept and each slope (in [-1, 1]).
#' @param slope_corr Numeric scalar or numeric matrix. If scalar, used as a common
#'   pairwise correlation among slopes; if matrix, must be `n_random x n_random`
#'   and provides slope–slope correlations.
#' @param sigma Numeric. Residual (observation-level) standard deviation.
#' @param seed Integer or NULL. Optional RNG seed for reproducibility.
#' @param random_vars Integer vector or NULL. Indices of the predictors that receive
#'   random slopes; if NULL random indices are selected.
#' @param intercept Numeric. Fixed intercept added to the linear predictor.
#'
#' @details
#' The function builds a `(1 + n_random) x (1 + n_random)` random-effects
#' covariance matrix where the first row/column correspond to the intercept.
#' Off-diagonal entries between intercept and slopes are set as
#' `intercept_slope_cor * intercept_sd * slope_sd_j`. Slopes are uncorrelated
#' by default unless `slope_corr` is supplied. The covariance matrix is checked
#' for positive-definiteness and an error is raised if it is not valid.
#'
#' @return A data.frame with columns: \code{id}, \code{group} (factor), \code{y},
#' and simulated predictors \code{X1}..\code{Xp}.
#'
#' @examples
#' 
#' df <- generate_mixed_sim(
#'   n = 500, n_groups = 25, n_vars = 8, n_random = 2,
#'   rho = 0.4, random_slope_sd = c(0.4, 0.6), intercept_sd = 0.8,
#'   intercept_slope_cor = c(0.1, -0.2), slope_corr = 0.1, seed = 42
#' )
#' head(df)
#' 
#'
#' @importFrom MASS mvrnorm
#' @export
generate_mixed_sim <- function(
  n = 1000,
  n_groups = 20,
  n_vars = 10,
  n_random = 3,
  rho = 0.3,
  beta = NULL,
  mu = 0,                       # scalar or numeric vector length = n_vars
  random_slope_sd = 0.5,           # scalar or numeric vector length = n_random
  intercept_sd = 0.8,
  intercept_slope_cor = 0.1,       # scalar or numeric vector length = n_random, values in [-1,1]
  slope_corr = 0.1,                # optional: scalar (pairwise) or matrix for slope-slope correlations (default 0)
  sigma = 1.0,
  seed = NULL,
  random_vars = NULL,
  intercept = 0                     # fixed intercept
) {
  if (!is.null(seed)) set.seed(seed)
  if (n_random > n_vars) stop("n_random must be <= n_vars")
  if (rho < 0 || rho >= 1) stop("rho must be in [0, 1)")
  if (!is.null(random_vars)) {
    if (!all(random_vars %in% seq_len(n_vars))) stop("random_vars must be valid indices")
    n_random <- length(random_vars)
  }
  if (length(mu) == 1) mu <- rep(mu, n_vars)
  if (length(mu) != n_vars) stop("mu must be NULL, a scalar, or a numeric vector of length n_vars")
  # validate random_slope_sd
  if (length(random_slope_sd) == 1) random_slope_sd <- rep(random_slope_sd, n_random)
  if (length(random_slope_sd) != n_random) stop("random_slope_sd must be scalar or length n_random")
  if (length(intercept_slope_cor) == 1) intercept_slope_cor <- rep(intercept_slope_cor, n_random)
  if (length(intercept_slope_cor) != n_random) stop("intercept_slope_cor must be scalar or length n_random")
  if (any(intercept_slope_cor < -1 | intercept_slope_cor > 1)) stop("intercept_slope_cor must be in [-1, 1]")
  # group assignment (balanced as possible)
  groups <- rep(seq_len(n_groups), length.out = n)
  groups <- sample(groups, size = n, replace = FALSE)
  # predictor covariance (compound symmetry)
  Sigma <- matrix(rho, n_vars, n_vars)
  diag(Sigma) <- 1
  Xmat <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  colnames(Xmat) <- paste0("X", seq_len(n_vars))
  # fixed-effect coefficients
  if (is.null(beta)) {
    beta <- stats::rnorm(n_vars, mean = 0, sd = 1)
  } else {
    if (length(beta) != n_vars) stop("length(beta) must equal n_vars")
  }
  # choose which predictors have random slopes
  if (is.null(random_vars)) {
    random_vars <- sample(seq_len(n_vars), size = n_random)
  } else {
    n_random <- length(random_vars)
  }
  # construct random-effects covariance (1 + n_random) x (1 + n_random)
  re_dim <- 1 + n_random
  Sigma_re <- matrix(0, nrow = re_dim, ncol = re_dim)
  # variances
  Sigma_re[1, 1] <- intercept_sd^2
  Sigma_re[2:re_dim, 2:re_dim] <- diag(random_slope_sd^2)
  # intercept <-> slope covariances
  Sigma_re[1, 2:re_dim] <- intercept_slope_cor * intercept_sd * random_slope_sd
  Sigma_re[2:re_dim, 1] <- Sigma_re[1, 2:re_dim]
  # slopes pairwise correlations (optional)
  if (length(slope_corr) == 1) {
    if (slope_corr != 0) {
      off <- slope_corr
      for (i in 2:re_dim) for (j in 2:re_dim) if (i != j) {
        Sigma_re[i, j] <- off * sqrt(Sigma_re[i, i] * Sigma_re[j, j])
      }
    }
  } else {
    slope_corr <- as.matrix(slope_corr)
    if (!all(dim(slope_corr) == c(n_random, n_random))) stop("slope_corr must be scalar or n_random x n_random matrix")
    for (i in seq_len(n_random)) for (j in seq_len(n_random)) {
      Sigma_re[i + 1, j + 1] <- slope_corr[i, j] * random_slope_sd[i] * random_slope_sd[j]
    }
  }
  # check PD
  ev <- eigen(Sigma_re, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 1e-8) stop("constructed random-effects covariance is not positive-definite; adjust SDs/correlations")
  # draw group-level random effects (rows = groups, cols = intercept + slopes)
  re_draws <- MASS::mvrnorm(n = n_groups, mu = rep(0, re_dim), Sigma = Sigma_re)
  colnames(re_draws) <- c("(Intercept)", paste0("bX", random_vars))
  rownames(re_draws) <- paste0("G", seq_len(n_groups))
  # build response
  fixed_part <- as.vector(Xmat %*% beta) + intercept
  rand_contrib <- numeric(n)
  intercept_part <- numeric(n)
  for (i in seq_len(n)) {
    g <- groups[i]
    intercept_part[i] <- re_draws[g, 1]
    xis <- Xmat[i, random_vars, drop = FALSE]
    rand_coefs <- re_draws[g, 2:re_dim, drop = FALSE]
    rand_contrib[i] <- sum(xis * rand_coefs)
  }
  eps <- stats::rnorm(n, mean = 0, sd = sigma)
  y <-  fixed_part + intercept_part + rand_contrib + eps
  df <- data.frame(
    id = seq_len(n),
    group = factor(paste0("G", groups)),
    y = y,
    Xmat,
    stringsAsFactors = FALSE
  )
  return(df)
}