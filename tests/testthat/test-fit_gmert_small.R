test_that("fit_gmert_small returns the expected object structure (fast settings)", {
  # NOTE: fit_gmert_small is computationally expensive (nested EM/PQL loops),
  # even if optimized via A^{-1}. We cap iterations to keep tests fast.
  # We cap both outer and inner iterations to keep the test fast.
  # This shouldn't be done in practice! (the GLL may not converge properly.)
  df <- sim_data_gmert(G = 10, n_i = 25, seed = 1)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 1)

  fit <- fit_gmert_small(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 20,
    max_iter_out = 20,
    tol = 1e-4,
    maxdepth = 3,
    minsplit = 20,
    minbucket = 10,
    xval = 0,
    save_tree = TRUE,
    save_gll = TRUE,
    save_train_ids = TRUE,
    sanity_checks = FALSE
  )

  # basic structure
  expect_type(fit, "list")
  expect_true(all(c(
    "tree", "b", "D", "sigma2", "mu",
    "converged_in", "converged_out",
    "n_iter", "train_ids", "gll", "tol" ) %in% names(fit)))

  # fitted tree
  expect_true(inherits(fit$tree, "rpart"))

  # random effects matrix dimensions
  expect_true(is.matrix(fit$b))
  expect_equal(ncol(fit$b), 2) # intercept + slope on x1
  expect_equal(nrow(fit$b), length(unique(split$train$id)))

  # covariance matrix D
  expect_true(is.matrix(fit$D))
  expect_equal(dim(fit$D), c(2, 2))

  # sigma2, tol, time
  expect_true(is.numeric(fit$sigma2))
  expect_length(fit$sigma2, 1)
  expect_true(is.finite(fit$sigma2))

  expect_true(is.numeric(fit$tol))
  expect_length(fit$tol, 1)


  # mu: probabilities for training rows
  expect_true(is.numeric(fit$mu))
  expect_equal(length(fit$mu), nrow(split$train))
  expect_true(all(is.finite(fit$mu)))
  expect_true(all(fit$mu >= 0 & fit$mu <= 1))

  # convergence flags
  expect_true(is.logical(fit$converged_out))
  expect_length(fit$converged_out, 1)
  expect_true(is.logical(fit$converged_in))

  # stored train ids
  expect_equal(length(fit$train_ids), nrow(split$train))

  # gll trace exists
  expect_true(is.numeric(fit$gll))
  expect_true(length(fit$gll) >= 2)
})

test_that("fit_gmert_small runs end-to-end and predicts without error", {
  # NOTE: capped iterations for speed in tests.
  df <- sim_data_gmert(G = 8, n_i = 20, seed = 2)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 2)

  fit <- fit_gmert_small(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 15,
    max_iter_out = 15,
    tol = 1e-4,
    maxdepth = 3,
    minsplit = 15,
    minbucket = 8,
    xval = 0,
    save_tree = TRUE,
    save_gll = TRUE,
    save_train_ids = TRUE,
    sanity_checks = FALSE
  )

  pred <- predict_gmert(
    fit = fit,
    new_df = split$test,
    random_effect = "x1",
    thr = 0.5,
    id = "id"
  )

  expect_true(is.numeric(pred) || is.integer(pred))
  expect_equal(length(pred), nrow(split$test))
  expect_true(all(pred %in% c(0, 1)))
})

test_that("fit_gmert_small is robust to different seeds/data (smoke test)", {
  # NOTE: capped iterations for speed in tests.
  df <- sim_data_gmert(G = 6, n_i = 15, seed = 99)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 99)

  expect_error(
    fit_gmert_small(
      df = split$train,
      id = "id",
      target = "y",
      random_effects = "x1",
      max_iter_inn = 10,
      max_iter_out = 10,
      tol = 1e-4,
      maxdepth = 2,
      minsplit = 10,
      minbucket = 5,
      xval = 0
    ),
    NA
  )
})

test_that("fit_gmert_small and fit_gmert converge to similar results (same data/seed)", {
  # NOTE: These functions are computationally expensive (nested EM/PQL loops).
  # We keep max_iter_out <= 20 for test runtime, but allow many inner iterations
  # so both methods can get close to the same fixed point.
  df <- sim_data_gmert(G = 8, n_i = 20, seed = 123)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 123)

  # Try to reduce any randomness in the tree fit:
  # - set.seed before each fit
  # - set xval = 0 (no CV)
  # - keep tree hyperparams identical
  set.seed(999)
  fit_std <- fit_gmert(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 500,
    max_iter_out = 20,
    tol = 1e-6,
    maxdepth = 3,
    minsplit = 15,
    minbucket = 8,
    xval = 0
  )

  set.seed(999)
  fit_small <- fit_gmert_small(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 500,
    max_iter_out = 20,
    tol = 1e-6,
    maxdepth = 3,
    minsplit = 15,
    minbucket = 8,
    xval = 0,
    save_tree = TRUE
  )

  # --- Similarity checks (tolerances are intentionally not ultra-tight) ---
  # mu: mean absolute difference should be small
  mad_mu <- mean(abs(fit_std$mu - fit_small$mu))
  expect_true(is.finite(mad_mu))
  expect_lt(mad_mu, 0.05)

  # sigma2 should be close
  expect_true(is.finite(fit_std$sigma2))
  expect_true(is.finite(fit_small$sigma2))
  expect_lt(abs(fit_std$sigma2 - fit_small$sigma2), 0.1)

  # D should be close (Frobenius norm)
  D_diff <- sqrt(sum((fit_std$D - fit_small$D)^2))
  expect_true(is.finite(D_diff))
  expect_lt(D_diff, 0.2)

  # random effects b should be close on average
  b_mad <- mean(abs(fit_std$b - fit_small$b))
  expect_true(is.finite(b_mad))
  expect_lt(b_mad, 0.2)

  # predicted classes on training should mostly agree
  yhat_std   <- predict_gmert(fit_std,   new_df = split$train, random_effect = "x1", thr = 0.5, id = "id")
  yhat_small <- predict_gmert(fit_small, new_df = split$train, random_effect = "x1", thr = 0.5, id = "id")
  agree <- mean(yhat_std == yhat_small)
  expect_true(is.finite(agree))
  expect_gt(agree, 0.9)
})

