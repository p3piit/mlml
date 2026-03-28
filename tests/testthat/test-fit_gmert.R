test_that("fit_gmert returns the expected object structure (fast settings)", {
  # NOTE: fit_gmert is computationally expensive (nested EM/PQL loops).
  # We cap both outer and inner iterations to keep the test fast.
  # This shouldn't be done in practice! (the GLL may not converge properly.)
  df <- sim_data_gmert(G = 10, n_i = 25, seed = 1)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 1)

  fit <- fit_gmert(
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
    xval = 0
  )

  # basic structure
  expect_type(fit, "list")
  expect_true(all(c(
    "tree", "b", "D", "sigma2", "mu",
    "converged_in", "converged_out",
    "n_iter", "train_ids", "gll", "tol"
  ) %in% names(fit)))

  # fitted tree
  expect_true(inherits(fit$tree, "rpart"))

  # random effects matrix dimensions
  expect_true(is.matrix(fit$b))
  expect_equal(ncol(fit$b), 2) # intercept + slope on x1
  expect_equal(nrow(fit$b), length(unique(split$train$id)))

  # covariance matrix D
  expect_true(is.matrix(fit$D))
  expect_equal(dim(fit$D), c(2, 2))

  # sigma2 and tol
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

test_that("fit_gmert runs end-to-end and predicts without error", {
  # NOTE: capped iterations for speed in tests.
  df <- sim_data_gmert(G = 8, n_i = 20, seed = 2)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 2)

  fit <- fit_gmert(
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
    xval = 0
  )

  # Predictions on the test set
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

test_that("fit_gmert is robust to different seeds/data (smoke test)", {
  # NOTE: capped iterations for speed in tests.
  df <- sim_data_gmert(G = 6, n_i = 15, seed = 99)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 99)

  expect_error(
    fit_gmert(
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

