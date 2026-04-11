# Main fitting function for categorical outcomes using Ainv-based updates

#' Fit a generalized mixed-effects regression tree (GMERT) using A-inverse updates
#'
#' Extends the binary GMERT to K-class categorical responses using the multinomial
#' logit link. The fixed part is represented by K-1 regression trees (one per
#' log-odds contrast against the reference class K), while the random part retains
#' a linear mixed-effects structure with a covariance matrix of dimension
#' (K-1)*q x (K-1)*q.
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
#' @param target Character scalar. Name of the categorical response variable in
#'   \code{df}. Must be a factor or integer with values in \{1, ..., K\}.
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
#' @param sanity_checks Logical. Whether to print sanity check messages during fitting.
#' @param initial_mu Optional initial class probabilities. If NULL, uniform 1/K is used.
#'
#' @return A list with fitted model components.
#'
#' @export
#'
fit_gmert_cat <- function(df,
                            id = "id",
                            target = "y",
                            random_effects = "x1",
                            max_iter_inn = 1000,
                            max_iter_out = 1000,
                            tol = 1e-6,
                            cp = 0.0,
                            minsplit = 50,
                            minbucket = 20,
                            maxdepth = 5,
                            xval = 10,
                            sanity_checks = FALSE,
                            initial_mu = NULL) {

  # --- Basic setup ---
  N <- nrow(df)
  G <- length(unique(df[[id]]))
  idx_by_cluster <- split(seq_len(N), df[[id]]) %>%
    .[lengths(.) > 0]

  y_raw <- df[[target]]

  # Encode response as integer class labels 1..K
  if (is.factor(y_raw)) {
    classes <- levels(y_raw)
    y_int <- as.integer(y_raw)
  } else {
    classes <- as.character(sort(unique(y_raw)))
    y_int <- match(as.character(y_raw), classes)
  }

  K  <- length(classes)
  K1 <- K - 1L

  if (K < 3L) {
    stop("fit_gmert_small requires K >= 3. For binary outcomes use the binary version.")
  }

  # One-hot encoding: N x K
  Y_oh <- matrix(0, N, K)
  Y_oh[cbind(seq_len(N), y_int)] <- 1

  # Random-effects design matrix
  Z <- as.matrix(cbind(1, df[random_effects]))
  q <- ncol(Z)

  # Fixed-effects data used by rpart
  Xdf <- data.table::as.data.table(
    df[setdiff(names(df), c(id, target))]
  )

  # Tree controls
  ctrl <- rpart.control(cp = cp,
                        minsplit = minsplit,
                        xval = xval,
                        minbucket = minbucket,
                        maxdepth = maxdepth)

  # --- Initialization (Step 0) ---
  M <- 0L

  # Initial class probabilities: N x K
  if (is.null(initial_mu)) {
    mu <- matrix(1 / K, nrow = N, ncol = K)
  } else if (length(initial_mu) == 1L) {
    mu <- matrix(initial_mu, nrow = N, ncol = K)
    mu <- mu / rowSums(mu)
  } else if (length(initial_mu) == K) {
    mu <- matrix(rep(initial_mu, each = N), nrow = N, ncol = K)
    mu <- mu / rowSums(mu)
  } else if (is.matrix(initial_mu) && all(dim(initial_mu) == c(N, K))) {
    mu <- initial_mu / rowSums(initial_mu)
  } else {
    stop("initial_mu must be NULL, a scalar, a length-K vector, or an N x K matrix.")
  }

  # Numerical safeguard
  mu <- pmax(mu, 1e-15)
  mu <- mu / rowSums(mu)

  # Initial linear predictor for non-reference classes: N x (K-1)
  eta <- log(mu[, 1:K1, drop = FALSE] / mu[, K])

  # Initial pseudo-response and working weights
  y_t <- pseudo_fun_cat(Y = Y_oh,
                        mu = mu,
                        eta = eta,
                        N = N,
                        K1 = K1)

  W <- W_fun_cat(mu = mu,
                 idx = idx_by_cluster)

  sigma2 <- 1
  D <- diag(K1 * q)
  b <- matrix(0, G, K1 * q)

  gll <- c(0)
  eta_old <- matrix(0, N, K1)
  converged_in <- c()
  converged_out <- FALSE
  d_eta <- numeric(max_iter_out)
  time_start <- proc.time()

  # --- Outer loop (PQL updates) ---
  repeat {
    m = 0                                       # reset inner-loop counter

    # --- Inner loop (EM-like iterations) ---
        # --- Inner loop (EM-like iterations) ---
    repeat {
      m <- m + 1L

      # (1.i) Partial E-step: compute adjusted pseudo-response
      # y_star = y_t - Z_tilde b, stored as N x K1
      zb <- matrix(0, N, K1)

      for (g in seq_len(G)) {
        idg <- idx_by_cluster[[g]]
        Zi  <- Z[idg, , drop = FALSE]

        for (k in seq_len(K1)) {
          cols_k <- ((k - 1L) * q + 1L):(k * q)
          zb[idg, k] <- Zi %*% b[g, cols_k]
        }
      }

      y_star <- y_t - zb   # N x K1

      # (1.ii) M-step: fit one regression tree per contrast
      trees <- vector("list", K1)
      fhat  <- matrix(0, N, K1)

      for (k in seq_len(K1)) {
        Xdf[, y_star_k := y_star[, k]]

        trees[[k]] <- rpart(y_star_k ~ .,
                            data = Xdf,
                            method = "anova",
                            control = ctrl)

        fhat[, k] <- as.numeric(predict(trees[[k]], newdata = Xdf))

        Xdf[, y_star_k := NULL]
      }

      # (1.iii) Update random effects b_i
      Ainv <- Ajnv_fun_cat(G = G,
                           Z = Z,
                           W = W,
                           D = D,
                           sigma2 = sigma2,
                           idx = idx_by_cluster,
                           K1 = K1)

      b <- b_fun_cat(G = G,
                           Z = Z,
                           W = W,
                           idx = idx_by_cluster,
                           y_t = y_t,
                           fhat = fhat,
                           Ainv = Ainv,
                           sigma2 = sigma2,
                           K1 = K1)

      # (2.i) Update sigma^2
      sigma2 <- sigma_fun_cat(N = N,
                                    G = G,
                                    idx = idx_by_cluster,
                                    Z = Z,
                                    W = W,
                                    y_t = y_t,
                                    fhat = fhat,
                                    b = b,
                                    Ainv = Ainv,
                                    K1 = K1)

      # (2.ii) Update D
      D <- D_fun_cat(G = G,
                           b = b,
                           Ainv = Ainv)

      # --- Inner-loop convergence check (GLL stabilization) ---
      gll[m + 1] <- gll_fun_cat(idx = idx_by_cluster,
                                b = b,
                                y = y_t,
                                Z = Z,
                                D = D,
                                s2 = sigma2,
                                fhat = fhat,
                                W = W,
                                K1 = K1)

      if (m > 1L) {
        rel <- abs(gll[m + 1] - gll[m]) / (abs(gll[m]) + 1e-12)
        if (rel < tol) {
          n_iter <- m
          converged_in_t <- TRUE
          break
        }
      }

      if (m >= max_iter_inn) {
        n_iter <- m
        converged_in_t <- FALSE
        break
      }
    }

    # --- Outer-loop update (PQL step) ---
    converged_in <- c(converged_in, converged_in_t)

    zb <- matrix(0, N, K1)
    for (g in seq_len(G)) {
      idg <- idx_by_cluster[[g]]
      Zi  <- Z[idg, , drop = FALSE]

      for (k in seq_len(K1)) {
        cols_k <- ((k - 1L) * q + 1L):(k * q)
        zb[idg, k] <- Zi %*% b[g, cols_k]
      }
    }

    eta <- fhat + zb   # N x K1

    # Convert eta to class probabilities mu via softmax with class K as reference
    exp_eta <- exp(eta)
    exp_eta[is.infinite(exp_eta)] <- 1e50
    denom <- 1 + rowSums(exp_eta)

    mu <- matrix(0, N, K)
    mu[, 1:K1] <- exp_eta / denom
    mu[, K] <- 1 / denom

    # Numerical safeguard
    mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
    mu <- mu / rowSums(mu)

    M <- M + 1L

    # Outer stopping rule
    d_eta[M] <- sqrt(mean((eta - eta_old)^2))
    if (d_eta[M] < tol) {
      converged_out <- TRUE
      break
    }

    if (M >= max_iter_out) {
      converged_out <- FALSE
      message(sprintf("WARNING: the PQL algorithm did not converge in %d iterations.",
                      max_iter_out))
      break
    }

    # Update working quantities for next outer iteration
    eta_old <- eta

    y_t <- pseudo_fun_cat(Y = Y_oh,
                          mu = mu,
                          eta = eta,
                          N = N,
                          K1 = K1)

    W <- W_fun_cat(mu = mu,
                   idx = idx_by_cluster)

    if (M %% 10 == 0 && sanity_checks) {
      time_elapsed <- proc.time() - time_start
      time_elapsed_seconds <- time_elapsed["elapsed"]
      time_elapsed_minutes <- floor(time_elapsed["elapsed"] / 60) -
        floor(time_elapsed["elapsed"] / 3600) * 60
      time_elapsed_hours <- floor(time_elapsed["elapsed"] / 3600)

      print(sprintf(
        paste0("Outer iteration %d: elapsed time = %.0f hours, %.0f minutes, %.2f seconds,\n",
               "d_eta = %.6f.\n"),
        M,
        time_elapsed_hours,
        time_elapsed_minutes,
        time_elapsed_seconds - floor(time_elapsed_seconds / 60) * 60,
        d_eta[M]
      ))
    }
  }

  time_elapsed <- proc.time() - time_start
  time_elapsed_seconds <- time_elapsed["elapsed"] -
    floor(time_elapsed["elapsed"] / 60) * 60
  time_elapsed_minutes <- floor(time_elapsed["elapsed"] / 60) -
    floor(time_elapsed["elapsed"] / 3600) * 60
  time_elapsed_hours <- floor(time_elapsed["elapsed"] / 3600)

  message(sprintf("Total elapsed time: %.0f hours, %.0f minutes, %.2f seconds.\n",
                  time_elapsed_hours, time_elapsed_minutes, time_elapsed_seconds))

  out <- list(
    trees = trees,
    b = b,
    D = D,
    sigma2 = sigma2,
    mu = mu,
    K = K,
    classes = classes,
    converged_in = converged_in,
    converged_out = converged_out,
    n_iter = n_iter,
    tol = tol,
    train_ids = unique(df[[id]]),
    d_eta = d_eta[seq_len(M)]
  )

  rm(Xdf)
  gc()

  return(out)
}
