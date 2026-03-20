test_that("fit_gmerf returns a valid fitted object and makes predictions", {
  # NOTE: fit_gmerf is computationally expensive (nested outer/inner loops + random forest).
  # We cap both outer and inner iterations to keep the test fast.
  # This shouldn't be done in practice! (the GLL may not converge properly.)
  df <- sim_data_gmert(G = 6, n_i = 20, seed = 42)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 42)

  set.seed(123)
  fit <- fit_gmerf(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 200,
    max_iter_out = 20,
    tol = 1e-4,
    ntrees = 50,         # keep small for tests
    min_node_size = 5,
    mtry = NULL,
    max.depth = NULL,
    seed = 123
  )

  # structure checks
  expect_type(fit, "list")
  expect_true(all(c("forest", "b", "D", "sigma2", "mu", "converged_in",
                    "converged_out", "n_iter", "train_ids", "gll", "tol") %in% names(fit)))

  # sanity checks on returned objects
  expect_true(is.list(fit$forest) || is.environment(fit$forest) || !is.null(fit$forest))
  expect_true(is.matrix(fit$b))
  expect_true(is.matrix(fit$D))
  expect_true(length(fit$mu) == nrow(split$train))
  expect_true(is.numeric(fit$sigma2) && length(fit$sigma2) == 1 && is.finite(fit$sigma2) && fit$sigma2 > 0)

  # dimensions
  G_train <- length(unique(split$train[[ "id" ]]))
  expect_equal(nrow(fit$b), G_train)
  expect_equal(ncol(fit$D), nrow(fit$D))
  expect_equal(ncol(fit$b), nrow(fit$D))

  # predictions
  yhat <- predict_gmerf(fit, new_df = split$test, thr = 0.5, random_effects = "x1", id = "id")
  expect_length(yhat, nrow(split$test))
  expect_true(all(yhat %in% c(0, 1)))
})

test_that("fit_gmerf is reproducible given the same seed", {
  # NOTE: fit_gmerf includes randomness from the random forest (and possibly upstream).
  # We use a small number of trees to keep the test fast.
  df <- sim_data_gmert(G = 6, n_i = 20, seed = 7)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 7)

  set.seed(999)
  fit1 <- fit_gmerf(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 150,
    max_iter_out = 20,
    tol = 1e-4,
    ntrees = 50,
    seed = 999
  )

  set.seed(999)
  fit2 <- fit_gmerf(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 150,
    max_iter_out = 20,
    tol = 1e-4,
    ntrees = 50,
    seed = 999
  )

  # compare key numeric outputs (should be identical or extremely close)
  expect_equal(fit1$sigma2, fit2$sigma2, tolerance = 1e-10)
  expect_equal(fit1$D, fit2$D, tolerance = 1e-10)
  expect_equal(fit1$b, fit2$b, tolerance = 1e-10)
  expect_equal(fit1$mu, fit2$mu, tolerance = 1e-10)

  # predictions should match exactly
  yhat1 <- predict_gmerf(fit1, new_df = split$test, thr = 0.5, random_effects = "x1", id = "id")
  yhat2 <- predict_gmerf(fit2, new_df = split$test, thr = 0.5, random_effects = "x1", id = "id")
  expect_identical(yhat1, yhat2)
})
