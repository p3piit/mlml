#' Fit a (small/optimized) Generalized Mixed-Effects Random Forest model (GMERF)
#'
#' Optimized variant of fit_gmerf_small for large simulations:
#' - forces single-threaded BLAS / MKL / OpenMP inside the function
#' - precomputes the RF data frame once (avoids repeated subsetting)
#' - replaces per-cluster R-loops used to build Z %*% b with a vectorized mapping
#' - forces ranger to use a single thread (num.threads = 1)
#'
#' The original implementation lines that were replaced are kept as commented blocks
#' immediately above each change to ease comparison.
#'
#' @inheritParams fit_gmerf_small
#' @export
fit_gmerf_opt <- function(df,               # df: data.frame with columns
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
                            seed = 1234,          # random seed for reproducibility
                            num.threads = 1     # number of threads for ranger
) {

  # --- Force single-threading for BLAS / OpenMP / MKL inside this function ---
  old_openblas <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  old_omp      <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  old_mkl      <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

  on.exit({
    # restore previous env values (if any)
    if (!is.na(old_openblas)) Sys.setenv(OPENBLAS_NUM_THREADS = old_openblas) else Sys.unsetenv("OPENBLAS_NUM_THREADS")
    if (!is.na(old_omp))      Sys.setenv(OMP_NUM_THREADS = old_omp)      else Sys.unsetenv("OMP_NUM_THREADS")
    if (!is.na(old_mkl))      Sys.setenv(MKL_NUM_THREADS = old_mkl)      else Sys.unsetenv("MKL_NUM_THREADS")
  }, add = TRUE)

  # --- Basic setup ---
  N <- nrow(df)                                 # total number of observations
  idx_by_cluster <- split(seq_len(N), df[[id]])    # list: row indices grouped by cluster
  G <- length(idx_by_cluster)                    # number of clusters

  y <- df[[target]]                                     # response vector
  Z <- cbind(1, df[[random_effects]])                          # random-effects design: intercept + slope on x1
  q <- ncol(Z)                                  # number of random effects (q = 2)

  # Precompute the RF covariate set once (avoid repeated subsetting inside inner loop).
  # OLD:
  #  # (original did this inside the inner loop)
  #  # Xrf <- df[setdiff(names(df), c(id, target))]
  Xrf_cols <- setdiff(names(df), c(id, target))
  Xrf <- df[Xrf_cols]

  # Precompute a vector mapping each observation to its cluster index (matching idx_by_cluster order).
  # This will be used to vectorize Z %*% b evaluations by indexing b by-observation.
  # OLD:
  #  # zb <- numeric(N)
  #  # for (g in seq_len(G)) {
  #  #   idx <- idx_by_cluster[[g]]
  #  #   zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]
  #  # }
  cluster_names <- names(idx_by_cluster)
  cluster_index <- match(df[[id]], cluster_names)  # integer vector length N mapping each row to 1..G

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

  # --- Outer loop (PQL updates) ---
  repeat {
    m = 0                                       # reset inner-loop counter

    # --- Inner loop (EM-like iterations) ---
    repeat {
      m <- m + 1L

      # (1.i) Partial E-step: compute adjusted pseudo-response y_tilde* = y_tilde - Z b
      # OLD:
      #   zb <- numeric(N)
      #   for (g in seq_len(G)) {
      #     idx <- idx_by_cluster[[g]]              # indices for cluster g
      #     zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]  # Z_i b_i
      #   }
      # Vectorized: build b for each observation by indexing b with cluster_index,
      # then compute row-wise inner product with Z.
      b_by_obs <- b[cluster_index, , drop = FALSE]   # N x q
      zb <- rowSums(Z * b_by_obs)                    # fast vectorized Z b per obs

      y_star <- y_t - zb                        # adjusted response for tree fit

      # (1.ii) M-step: fit random forest to (y_tilde*, X, W)
      # OLD:
      #  rf <- ranger::ranger(
      #    formula         = y_star ~ .,
      #    data            = cbind(y_star = y_star, Xrf),
      #    case.weights    = w,
      #    num.trees       = ntrees,
      #    min.node.size   = min_node_size,
      #    max.depth       = max.depth,
      #    respect.unordered.factors = TRUE,
      #    importance      = "none",
      #    write.forest    = TRUE,
      #    seed            = seed
      #  )
      # Use a single thread for ranger (num.threads = 1) to avoid interfering with parallel tasks.
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

    # OLD:
    #  zb <- numeric(N)
    #  for (g in seq_len(G)) {
    #    idx <- idx_by_cluster[[g]]
    #    zb[idx] <- Z[idx, , drop = FALSE] %*% b[g, ]
    #  }
    # Reuse the vectorized approach to recompute zb for the updated b.
    b_by_obs <- b[cluster_index, , drop = FALSE]
    zb <- rowSums(Z * b_by_obs)

    eta <- fhat + zb                            # recompute linear predictor

    # (Outer stopping rule – paper style)
    d_eta <- sqrt(mean((eta - eta_old)^2))      # RMS change of eta
    if (d_eta < tol) {
      converged_out <- TRUE
      break
    }

    # Update working quantities for next outer iteration
    eta_old <- eta
    mu <- exp(eta) / (1 + exp(eta))             # updated conditional means
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
    train_ids = df[[id]],          # cluster identifiers in training set
    gll = gll,                     # GLL trace (for diagnostics)
    tol = tol                      # convergence tolerance used
  )
}

