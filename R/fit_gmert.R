#'
#' Fit a generalized mixed-effects regression tree (GMERT)
#'
#' Fits a generalized mixed-effects regression tree model for clustered binary
#' outcomes using a PQL outer loop and an EM-like inner loop. The fixed part is
#' estimated via an `rpart` regression tree on a working pseudo-response, and the
#' random part is modeled with cluster-specific random intercepts and slopes
#' (currently implemented as an intercept plus one random slope).
#'
#' The function returns a list containing the fitted tree, estimated random
#' effects, covariance parameters, and convergence diagnostics.
#'
#' @param df A data frame containing the response, cluster identifier, and predictors.
#' @param id Character string. Name of the cluster identifier column in `df`.
#' @param target Character string. Name of the binary response column in `df` (coded 0/1).
#' @param random_effects Character string. Name of the predictor used for the random slope.
#'   The random-effects design is built as an intercept plus this column.
#' @param max_iter_inn Integer. Maximum number of EM-like inner-loop iterations.
#' @param max_iter_out Integer. Maximum number of PQL outer-loop iterations.
#' @param tol Numeric scalar. Convergence tolerance used for both inner and outer loops.
#' @param cp Numeric scalar. `rpart` complexity parameter (pruning threshold).
#' @param minsplit Integer. Minimum number of observations required to attempt a split in `rpart`.
#' @param minbucket Integer. Minimum number of observations in any terminal node in `rpart`.
#' @param maxdepth Integer. Maximum depth of the fitted `rpart` tree.
#' @param xval Integer. Number of cross-validation folds used by `rpart`.
#'
#' @return A list with components:
#' \describe{
#'   \item{tree}{An `rpart` object: the fitted regression tree for the fixed effects.}
#'   \item{b}{Numeric matrix (G x q) of estimated random effects by cluster.}
#'   \item{D}{Numeric matrix (q x q): estimated random-effects covariance matrix.}
#'   \item{sigma2}{Numeric scalar: estimated residual variance on the working scale.}
#'   \item{mu}{Numeric vector length N: fitted conditional means at convergence.}
#'   \item{converged_in}{Logical vector: inner-loop convergence flag for each outer iteration.}
#'   \item{converged_out}{Logical scalar: outer-loop convergence flag.}
#'   \item{n_iter}{Integer: number of inner iterations performed in the last outer iteration.}
#'   \item{train_ids}{Vector of cluster identifiers from the training data.}
#'   \item{gll}{Numeric vector: generalized log-likelihood trace used for diagnostics.}
#'   \item{tol}{Numeric scalar: tolerance used.}
#' }
#'
#' @details
#' This implementation assumes a binary response and uses a logistic link via a
#' penalized quasi-likelihood (PQL) linearization. The random-effects structure
#' is currently an intercept plus one slope specified by `random_effects`.
#'
#' @seealso [predict_gmert()] for prediction on new data.
#'
#' @family gmert
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
#' split <- split_gmert_data(df, train_prop = 0.7, seed = 1)
#'
#' # Fit GMERT
#' fit <- fit_gmert(split$train, id = "id", target = "y", random_effects = "x1",
#'                  max_iter_inn = 50, max_iter_out = 50, tol = 1e-4,
#'                  maxdepth = 3, minsplit = 20, minbucket = 10)
#'
#' # Predict on test set
#' pred <- predict_gmert(fit, split$test, random_effect = "x1", id = "id")
#' table(truth = split$test$y, pred = pred)
#'
#' @references
#' Hajjem, A., Larocque, D., & Bellavance, F. (2017).
#' \emph{Generalized mixed effects regression trees}.
#' Statistics & Probability Letters, 126, 114--118.
#' \doi{10.1016/j.spl.2017.02.033}

fit_gmert    <- function(df,               # df: data.frame with columns
                         id = "id",      # cluster identifier column name
                         target = "y",        # response column name
                         random_effects = "x1",  # random effects design column names (intercept + slope on x1)
                         max_iter_inn = 1000,  # maximum number of EM iterations (inner loop)
                         max_iter_out = 1000,  # maximum number of PQL iterations (outer loop)
                         tol = 1e-6,           # convergence tolerance for both loops (Aitken or relative diff)
                         cp = 0.0,             # rpart complexity parameter (pruning threshold)
                         minsplit = 50,        # minimum number of obs required to attempt a split
                         minbucket = 20,       # minimum number of obs in any terminal node
                         maxdepth = 5,         # maximum tree depth
                         xval = 10             # number of cross-validation folds in rpart
) {

  # --- Basic setup ---
  N <- nrow(df)                                 # total number of observations
  G <- length(unique(df[[id]]))                    # number of clusters
  idx_by_cluster <- split(seq_len(N), df[[id]])    # list: row indices grouped by cluster

  y <- df[[target]]                                     # response vector
  Z <- cbind(1, df[[random_effects]])                          # random-effects design: intercept + slope on x1
  q <- ncol(Z)                                  # number of random effects (q = 2)

  # --- Initialization (Step 0) ---
  M <- 0L                                       # outer-loop counter
  mu <- ifelse(y == 1, 0.75, 0.25)              # initial conditional means
  y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # initial pseudo-response (PQL linearization)
  w <- mu * (1 - mu)                            # initial working weights
  sigma2 <- 1                                   # initial residual variance
  D <- diag(2)                                  # initial random-effects covariance (identity)
  b <- matrix(0, G, q)                          # initialize cluster random effects
  gll <- c(0, 0)                                # GLL storage (2 slots for Aitken acceleration)
  eta_old <- rep(0, N)                          # previous eta for outer-loop convergence check
  converged_in <- c()                           # inner-loop convergence flags (per outer iteration)
  converged_out <- FALSE                        # outer-loop convergence flag

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

      # (1.ii) M-step: fit regression tree for fixed effects f(X)
      ctrl <- rpart.control(cp = cp, minsplit = minsplit, xval = xval,
                            minbucket = minbucket, maxdepth = maxdepth)
      Xdf <- df[setdiff(names(df), c(id, target))]
      tree <- rpart(y_star ~ .,
                    data = cbind(y_star = y_star, Xdf),
                    weights = w, method = "anova", control = ctrl)
      fhat <- as.numeric(predict(tree, newdata = Xdf))

      # (1.iii) Update random effects b_i
      Vi <- Vi_fun(Z = Z, D = D, s2 = sigma2, G = G, idx = idx_by_cluster, w = w)  # build V_i per cluster
      b <- b_fun(G = G, Z = Z, Vi = Vi, D = D, idx = idx_by_cluster,
                 y_t = y_t, fhat = fhat, w = w)                                   # update b_i estimates

      # (2.i) Update sigma^2 (residual variance)
      sigma2 <- sigma_fun(N = N, G = G, idx = idx_by_cluster,
                          b = b, y = y_t, Z = Z, D = D,
                          Vi = Vi, s2 = sigma2, fhat = fhat, w = w)

      # (2.ii) Update D (random-effects covariance)
      D <- D_fun(G = G, idx = idx_by_cluster,
                 Z = Z, D = D, Vi = Vi, b = b, w = w)

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
    mu <- pmin(pmax(plogis(eta), 10e-15), 1 - 1e-15)            # updated conditional means (capped below 1)

    # (Outer stopping rule – paper style)
    d_eta <- sqrt(mean((eta - eta_old)^2))      # RMS change of eta
    if (d_eta < tol) {
      converged_out <- TRUE
      break
    }


    M <- M + 1L
    if (M >= max_iter_out) {                    # guard against outer non-convergence
      converged_out <- FALSE
      message(sprintf("WARNING: the PQL algorithm did not converge in %d iterations.", max_iter_out))
      break
    }

    # Update working quantities for next outer iteration
    eta_old <- eta
    y_t <- log(mu / (1 - mu)) + (y - mu) / (mu * (1 - mu))  # new pseudo-response
    w <- mu * (1 - mu)                          # new working weights

  }

  # --- Return fitted components ---
  list(
    tree = tree,                   # fitted rpart tree (fixed-effects function f(X))
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
