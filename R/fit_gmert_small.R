# Main fitting function using Ainv-based updates

#' Fit a generalized mixed-effects regression tree (GMERT) using A-inverse updates
#'
#' Fits a generalized mixed-effects regression tree model for clustered binary data,
#' using an optimized set of updates based on precomputing \eqn{A_i^{-1}} to speed up
#' estimation of the random effects covariance and variance components.
#'
#' This routine follows the penalized quasi-likelihood (PQL) outer loop and an EM-like
#' inner loop described for GMERT, while replacing direct \eqn{V_i^{-1}} computations
#' with cluster-wise precomputed \eqn{A_i^{-1}} matrices.
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
#' @param cp Numeric. Complexity parameter passed to \code{rpart.control()}.
#' @param minsplit Integer. Minimum number of observations required to attempt a split.
#' @param minbucket Integer. Minimum number of observations in any terminal node.
#' @param maxdepth Integer. Maximum depth of the fitted tree.
#' @param xval Integer. Number of cross-validation folds in \code{rpart}.
#' @param save_tree Logical. Whether to save the full \code{rpart} tree object in the output.
#' @param save_gll Logical. Whether to save the generalized log-likelihood trace.
#' @param save_train_ids Logical. Whether to save the training cluster IDs.
#' @param sanity_checks Logical. Whether to print sanity check messages during fitting.
#'
#' @return A list with components:
#' \describe{
#'   \item{tree}{Fitted \code{rpart} regression tree for the fixed-effects component.}
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
#'   \item{time}{Elapsed runtime in seconds.}
#' }
#'
#' @references
#' Hajjem, Ahlem, François Bellavance, and Denis Larocque. 2011. “Mixed Effects Regression Trees for Clustered Data.”
#' \emph{Statistics & Probability Letters} 81 (4): 451–59. \doi{10.1016/j.spl.2010.12.003}.
#'
#' Hajjem, Ahlem, Denis Larocque, and François Bellavance. 2017. “Generalized Mixed Effects Regression Trees.”
#' \emph{Statistics & Probability Letters} 126 (July): 114–18. \doi{10.1016/j.spl.2017.02.033}.
#'
#' @seealso \code{\link{fit_gmert}} for the baseline GMERT fitting routine using \eqn{V_i}.
#'
#' @import rpart
#' @import ranger
#' @import MASS
#' @import tidyverse
#' @import lme4
#' @import parallel
#' @import ggplot2
#' @import performance
#' @import data.table
#'
#' @export
#'
#' @examples
#' # Simulate clustered binary data
#' df <- sim_data_gmert(G = 10, n_i = 30, seed = 1)
#'
#' # Fit the optimized GMERT (A^{-1}-based updates)
#' fit <- fit_gmert_small(df, id = "id", target = "y", random_effects = "x1",
#'                        max_iter_inn = 50, max_iter_out = 50, tol = 1e-4)
#'
#' # Predicted classes on training data
#' pred <- predict_gmert(fit, new_df = df, random_effect = "x1", id = "id")
fit_gmert_small    <- function(df,               # df: data.frame with columns
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
                               xval = 10,            # number of cross-validation folds in rpart
                               save_tree = FALSE,    # whether to save the full tree object
                               save_gll = FALSE,       # whether to save the GLL trace
                               save_train_ids = FALSE, # whether to save the training cluster IDs
                               sanity_checks = FALSE   # whether to print sanity check messages during fitting
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
  gll <- c(0   )                                # GLL storage (2 slots for Aitken acceleration)
  eta_old <- rep(0, N)                          # previous eta for outer-loop convergence check
  converged_in <- c()                           # inner-loop convergence flags (per outer iteration)
  converged_out <- FALSE                        # outer-loop convergence flag
  Xdf <- data.table::as.data.table(
        df[setdiff(names(df), c(id, target))])
  
  time_start <- proc.time()                        # start timer

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
      ctrl <- rpart.control(cp = cp, 
                            minsplit = minsplit, 
                            xval = xval,
                            minbucket = minbucket, 
                            maxdepth = maxdepth)
      Xdf <- Xdf[, y_star := ..y_star]
      tree <- rpart(y_star ~ .,
                    data = Xdf,
                    weights = w, 
                    method = "anova", 
                    control = ctrl)
      fhat <- as.numeric(predict(tree, newdata = Xdf))
      
      if (!save_tree) {
        rm(tree)
      }

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
      gll[m + 1] <- gll_fun(D = D, b = b, idx = idx_by_cluster,
                            Z = Z, y = y_t, fhat = fhat, s2 = sigma2, w = w)
      if (m > 1L) {
        rel <- abs(gll[m + 1] - gll[m]) / (abs(gll[m]) + 1e-12)
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

    if (M %% 10 == 0 & sanity_checks) {
     time_elapsed <- proc.time() - time_start
     time_elapsed_seconds <- time_elapsed["elapsed"]
     time_elapsed_minutes <- floor(time_elapsed_seconds / 60) - floor(time_elapsed_seconds / 3600) * 60
     time_elapsed_hours <- floor(time_elapsed_seconds / 3600)
     message(sprintf("Outer iteration %d: elapsed time = %.2f hours, %.0f minutes, %.2f seconds,\n
                     d_eta = %.6f. \n",
                     M, time_elapsed_hours, time_elapsed_minutes,
                     time_elapsed_seconds - floor(time_elapsed_seconds / 60) * 60, d_eta))
    }


  }

  
  time_elapsed <- proc.time() - time_start
  time_elapsed_seconds <- time_elapsed["elapsed"] - floor(time_elapsed["elapsed"] / 60) * 60
  time_elapsed_minutes <- floor(time_elapsed_seconds / 60) - floor(time_elapsed_seconds / 3600) * 60
  time_elapsed_hours <- floor(time_elapsed_seconds / 3600)
  message(sprintf("Total elapsed time: %.2f hours, %.0f minutes, %.2f seconds. \n",
                  time_elapsed_hours, time_elapsed_minutes, time_elapsed_seconds))

  out <- list(
    # include fitted forest only if explicitly requested and rf exists
    b = b,
    D = D,
    sigma2 = sigma2,
    mu = mu,
    converged_in = converged_in,
    converged_out = converged_out,
    n_iter = n_iter,
    tol = tol
  )
  if (save_tree) {
    out$tree <- tree
  }

  if (save_gll) {
    out$gll <- gll
  }

  if (save_train_ids) {
    out$train_ids <- df[[id]]
  }

  rm(Xdf)
  gc()

  return(out)
}
