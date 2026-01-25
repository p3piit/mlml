#' Fit a generalized mixed-effects random forest (GMERF) using A-inverse updates
#'
#' Fits a generalized mixed-effects model where the fixed-effects component is learned
#' via a random forest and the random effects account for clustered/longitudinal
#' dependence. Estimation proceeds via a PQL outer loop and an EM-like inner loop,
#' using precomputed cluster-specific \eqn{A_i^{-1}} matrices to update random effects
#' and variance components efficiently.
#'
#' @family gmert
#'
#' @param df A data.frame containing the response, cluster identifier, random-effect
#'   covariate(s), and fixed-effect predictors.
#' @param id Character scalar. Name of the clustering/grouping variable in \code{df}.
#' @param target Character scalar. Name of the binary response variable in \code{df}.
#' @param random_effects Character scalar. Name(s) of the variable(s) used for the
#'   random-effects design (intercept is added internally).
#' @param max_iter_inn Integer. Maximum number of EM-like iterations (inner loop).
#' @param max_iter_out Integer. Maximum number of PQL iterations (outer loop).
#' @param tol Numeric. Convergence tolerance used for both inner and outer loops.
#' @param ntrees Integer. Number of trees for the random forest.
#' @param mtry Integer or \code{NULL}. Number of candidate features tried at each split
#'   (passed to \code{ranger::ranger}; \code{NULL} uses \code{ranger}'s default).
#' @param min_node_size Integer. Minimum terminal node size for the forest.
#' @param max.depth Integer or \code{NULL}. Maximum tree depth for the forest
#'   (\code{NULL} = unlimited).
#' @param seed Integer. Random seed for reproducibility.
#' @param num.threads Integer or \code{NULL}. Number of threads used by \code{ranger}
#'   (\code{NULL} uses all available threads).
#'
#' @return A list with components:
#' \describe{
#'   \item{forest}{Fitted random forest object (stored as \code{rf$forest}).}
#'   \item{b}{Matrix of estimated random effects (clusters x q).}
#'   \item{D}{Estimated random-effects covariance matrix.}
#'   \item{sigma2}{Estimated residual variance.}
#'   \item{mu}{Fitted conditional means.}
#'   \item{converged_in}{Logical vector of inner-loop convergence flags (per outer iteration).}
#'   \item{converged_out}{Logical flag for outer-loop convergence.}
#'   \item{n_iter}{Number of inner-loop iterations performed in the last outer iteration.}
#'   \item{train_ids}{Cluster identifiers used for training.}
#'   \item{gll}{Generalized log-likelihood trace used for diagnostics.}
#'   \item{tol}{Tolerance used for convergence checks.}
#' }
#'
#' @references
#' Hajjem, Ahlem, François Bellavance, and Denis Larocque. 2011.
#' “Mixed Effects Regression Trees for Clustered Data.” \emph{Statistics & Probability Letters} 81 (4): 451–59.
#' \doi{10.1016/j.spl.2010.12.003}.
#'
#' Hajjem, Ahlem, Denis Larocque, and François Bellavance. 2017.
#' “Generalized Mixed Effects Regression Trees.” \emph{Statistics & Probability Letters} 126 (July): 114–18.
#' \doi{10.1016/j.spl.2017.02.033}.
#'
#' Pellagatti, Massimo, Chiara Masci, Francesca Ieva, and Anna M. Paganoni. 2021.
#' “Generalized Mixed-Effects Random Forest: A Flexible Approach to Predict University Student Dropout.”
#' \emph{Statistical Analysis and Data Mining: The ASA Data Science Journal} 14 (3): 241–57.
#' \doi{10.1002/sam.11505}.
#'
#' @seealso \code{\link{fit_gmerf}} for the baseline GMERF fitting routine;
#'   \code{\link{fit_gmert_small}} for the tree-based \eqn{A^{-1}} optimized variant.
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
#'
#' @examples
#' # Simulate clustered binary data
#' df <- sim_data_gmert(G = 10, n_i = 30, seed = 1)
#'
#' # Fit the optimized GMERF (A^{-1}-based updates)
#' fit <- fit_gmerf_small(df, id = "id", target = "y", random_effects = "x1",
#'                        max_iter_inn = 50, max_iter_out = 50, tol = 1e-4,
#'                        ntrees = 100, seed = 1)
#'
#' # Predicted classes on training data
#' pred <- predict_gmerf(fit, new_df = df, random_effects = "x1", id = "id")
fit_gmerf_small    <- function(df,               # df: data.frame with columns
                               id = "id",      # cluster identifier column name
                               target = "y",        # response column name
                               random_effects = "x1",  # random effects design column names (intercept + slope on x1)
                               max_iter_inn = 1000,  # maximum number of EM iterations (inner loop)
                               max_iter_out = 1000,  # maximum number of PQL iterations (outer loop)
                               tol = 1e-6,           # convergence tolerance for both loops (Aitken or relative diff)
                               ntrees = 500,        # RF: number of trees
                               mtry   = NULL,       # RF: features tried at each split (default = floor(p/3) if NULL)
                               min_node_size = 5,   # RF: terminal node size
                               max.depth = NULL,    # RF: optional max depth (NULL = unlimited)
                               seed = 1234,         # random seed for reproducibility
                               num.threads = NULL   # RF: number of threads (NULL = all available)
) {

  # --- Basic setup ---
  N <- nrow(df)                                 # total number of observations
  G <- length(unique(df[[id]]))                    # number of clusters
  idx_by_cluster <- split(seq_len(N), df[[id]])    # list: row indices grouped by cluster

  y <- df[[target]]                                     # response vector
  Z <- as.matrix(cbind(1, df[random_effects]))                          # random-effects design: intercept + slope on x1
  q <- ncol(Z)                                  # number of random effects (q = 2)

  # --- Initialization (Step 0) ---
  M <- 0L                                       # outer-loop counter
  mu <- ifelse(y == 1, 0.75, 0.25)              # initial conditional means
  y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # initial pseudo-response (PQL linearization)
  w <- mu * (1 - mu)                            # initial working weights
  sigma2 <- 1                                   # initial residual variance
  D <- diag(q)                                  # initial random-effects covariance (identity)
  b <- matrix(0, G, q)                          # initialize cluster random effects
  gll <- c(0, 0)                                # GLL storage (2 slots for Aitken acceleration)
  eta_old <- rep(0, N)                          # previous eta for outer-loop convergence check
  converged_in <- c()                           # inner-loop convergence flags (per outer iteration)
  converged_out <- FALSE                        # outer-loop convergence flag
  time_start <- proc.time()                    # start timer

  # --- Outer loop (PQL updates) ---
  repeat {
    m = 0                                       # reset inner-loop counter

    # --- Inner loop (EM-like iterations) ---
    repeat {
      m <- m + 1L

      # (1.i) Partial E-step: compute adjusted pseudo-response y_tilde* = y_tilde - Z b
      zb <- numeric(N)                          # cluster-specific random contributions
      for (g in seq_len(G)) {
        idx <- idx_by_cluster[[g]]              # indices for cluster g
        zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]  # Z_i b_i
      }
      y_star <- y_t - zb                        # adjusted response for tree fit

      # (1.ii) M-step: fit random forest to (y_tilde*, X, W)
      Xrf <- df[setdiff(names(df), c(id, target))]
      rf <- ranger::ranger(
        formula         = y_star ~ .,
        data            = cbind(y_star = y_star, Xrf),
        case.weights    = w,
        num.trees       = ntrees,
        min.node.size   = min_node_size,
        max.depth       = max.depth,
        respect.unordered.factors = TRUE,
        importance      = "none",
        write.forest    = TRUE,
        seed            = seed,
        num.threads     = num.threads
      )
      fhat <- as.numeric(predict(rf, data = Xrf)$predictions)

      # (1.iii) Update random effects b_i
      Ainv <- Ajnv_fun(Z = Z, D = D, sigma2 = sigma2, G = G, idx = idx_by_cluster, w = w)  # build V_i per cluster
      b <- b_fun_small(G = G, Z = Z, idx = idx_by_cluster,
                       y_t = y_t, fhat = fhat, w = w, sigma2 = sigma2, Ainv = Ainv)                                   # update b_i estimates

      # (2.i) Update sigma^2 (residual variance)
      sigma2 <- sigma_fun_small(N = N, G = G, idx = idx_by_cluster,
                                b = b, y = y_t, Z = Z,
                                Ainv = Ainv, sigma2 = sigma2, fhat = fhat, w = w)

      # (2.ii) Update D (random-effects covariance)
      D <- D_fun_small(G = G, b = b, Ainv = Ainv)

      # --- Inner-loop convergence check (GLL stabilization) ---
      gll[m + 2] <- gll_fun(D = D, b = b, idx = idx_by_cluster,
                            Z = Z, y = y_t, fhat = fhat, s2 = sigma2, w = w)
      if (m > 1L) {
        rel <- abs(gll[m + 2] - gll[m + 1]) / (abs(gll[m + 1]) + 1e-12)
        if (rel < tol) { n_iter <- m; converged_in_t <- TRUE; break }
      }

      if (m >= max_iter_inn) {                  # max-iteration guard
        n_iter <- m
        converged_in_t <- FALSE
        # message(sprintf("WARNING: the EM algorithm did not converge in %d iterations.", max_iter_inn))
        break
      }
    }

    # --- Outer-loop update (PQL step) ---
    converged_in <- c(converged_in, converged_in_t)
    zb <- numeric(N)
    for (g in seq_len(G)) {
      idx <- idx_by_cluster[[g]]
      zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]
    }
    eta <- fhat + zb                            # recompute linear predictor

    # (Outer stopping rule – paper style)
    d_eta <- sqrt(mean((eta - eta_old)^2))      # RMS change of eta
    if (d_eta < tol) {
      converged_out <- TRUE
      break
    }

    # Update working quantities for next outer iteration
    eta_old <- eta
    mu <- pmin(pmax(plogis(eta), 10e-15), 1 - 1e-15)            # updated conditional means (capped below 1)
    y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # new pseudo-response
    w <- mu * (1 - mu)                          # new working weights

    M <- M + 1L
    if (M >= max_iter_out) {                    # guard against outer non-convergence
      converged_out <- FALSE
      message(sprintf("WARNING: the PQL algorithm did not converge in %d iterations.", max_iter_out))
      break
    }
  }

  # --- Return fitted components ---
  list(
    forest = rf$forest,           # fitted regression forest
    b = b,                         # estimated random effects (G x q)
    D = D,                         # estimated random-effects covariance
    sigma2 = sigma2,               # estimated residual variance
    mu = mu,                       # conditional means
    converged_in = converged_in,   # convergence flags for inner loops
    converged_out = converged_out, # convergence flag for outer loop
    n_iter = n_iter,               # iterations performed (inner loop)
    train_ids = df$id,             # cluster identifiers in training set
    gll = gll,                     # GLL trace (for diagnostics)
    tol = tol                      # convergence tolerance used
  )
}
