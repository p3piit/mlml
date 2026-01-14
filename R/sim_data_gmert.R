#'
#' Simulate clustered binary data for GMERT models
#'
#' Generates synthetic clustered data from a logistic mixed-effects model
#' with random intercepts and random slopes. The simulated dataset is suitable
#' for fitting and testing generalized mixed effects regression tree (GMERT)
#' models.
#'
#' @param G Integer. Number of clusters.
#' @param n_i Integer. Number of observations per cluster.
#' @param beta0 Numeric. Fixed intercept.
#' @param beta1 Numeric. Fixed effect for \code{x1}.
#' @param beta2 Numeric. Fixed effect for \code{x2}.
#' @param beta3 Numeric. Fixed effect for \code{x3}.
#' @param sigma_b0 Numeric. Standard deviation of the random intercept.
#' @param sigma_b1 Numeric. Standard deviation of the random slope for \code{x1}.
#' @param rho Numeric. Correlation between random intercept and random slope.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A data.frame with columns \code{id}, \code{y}, \code{x1}, \code{x2},
#'   and \code{x3}.
#'
#' @examples
#' data <- sim_data_gmert()
#' head(data)
#'
#' @family gmert
#' @seealso \code{\link{fit_gmert}}, \code{\link{predict_gmert}}
#'
#' @import rpart
#' @import ranger
#' @import MASS
#' @import tidyverse
#' @import lme4
#' @import parallel
#' @import ggplot2
#' @import performance
#'
#' @export
sim_data_gmert <- function(G = 50,          # number of clusters
                           n_i = 40,        # observations per cluster
                           beta0 = 0.5,     # fixed intercept
                           beta1 = 1.2,     # fixed effect for x1
                           beta2 = -0.8,    # fixed effect for x2
                           beta3 = 0.6,     # fixed effect for x3
                           sigma_b0 = 0.8,  # SD of random intercept
                           sigma_b1 = 0.5,  # SD of random slope for x1
                           rho = 0.2,       # correlation between b0 and b1
                           seed = 123) {    # random seed for reproducibility

  set.seed(seed)

  # covariance matrix of random effects (b0_i, b1_i)
  D <- matrix(c(sigma_b0^2, rho * sigma_b0 * sigma_b1,
                rho * sigma_b0 * sigma_b1, sigma_b1^2), 2, 2)

  # generate random effects for G clusters
  b <- MASS::mvrnorm(G, mu = c(0, 0), Sigma = D)

  # cluster identifiers
  id <- rep(1:G, each = n_i)

  # generate covariates: x1 (uniform), x2 (normal), x3 (binary)
  x1 <- runif(G * n_i, -2, 2)
  x2 <- rnorm(G * n_i, 0, 1)
  x3 <- rbinom(G * n_i, 1, 0.5)

  # linear predictor: fixed part + random part (b0_i + b1_i * x1)
  eta <- numeric(G * n_i)
  for (g in seq_len(G)) {
    idx <- which(id == g)
    eta[idx] <- beta0 + beta1 * x1[idx] + beta2 * x2[idx] + beta3 * x3[idx] +
      b[g, 1] + b[g, 2] * x1[idx]
  }

  # convert to probabilities via logistic link
  p <- 1 / (1 + exp(-eta))

  # binary outcome drawn from Bernoulli(p)
  y <- rbinom(G * n_i, 1, p)

  # return as data.frame ready for model fitting
  data.frame(
    id = factor(id),
    y = y,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
}
