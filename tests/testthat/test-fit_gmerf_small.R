test_that("fit_gmerf_small returns a valid fitted object and makes predictions", {
  # NOTE: fit_gmerf_small is computationally expensive (nested outer/inner loops + random forest),
  # even with the A^{-1}-based updates.
  # We cap both outer and inner iterations to keep the test fast.
  # This shouldn't be done in practice! (the GLL may not converge properly.)
  df <- sim_data_gmert(G = 6, n_i = 20, seed = 42)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 42)

  set.seed(123)
  fit <- fit_gmerf_small(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 200,
    max_iter_out = 20,
    tol = 1e-4,
    ntrees = 50,           # keep small for tests
    min_node_size = 5,
    mtry = NULL,
    max.depth = NULL,
    seed = 123,
    num.threads = 1,        # keep deterministic / light for checks
    sanity_checks = FALSE
  )

  # structure checks
  expect_type(fit, "list")
  expect_true(all(c("forest", "b", "D", "sigma2", "mu", "converged_in",
                    "converged_out", "n_iter", "train_ids", "tol") %in% names(fit)))

  # sanity checks on returned objects
  expect_true(is.list(fit$forest) || is.environment(fit$forest) || !is.null(fit$forest))
  expect_true(is.matrix(fit$b))
  expect_true(is.matrix(fit$D))
  expect_true(length(fit$mu) == nrow(split$train))
  expect_true(is.numeric(fit$sigma2) && length(fit$sigma2) == 1 && is.finite(fit$sigma2) && fit$sigma2 > 0)

  # dimensions
  G_train <- length(unique(split$train[["id"]]))
  expect_equal(nrow(fit$b), G_train)
  expect_equal(ncol(fit$D), nrow(fit$D))
  expect_equal(ncol(fit$b), nrow(fit$D))

  # predictions (use the exported predict_gmerf, since the "small" variant shares the same predictor)
  yhat <- predict_gmerf(fit, new_df = split$test, thr = 0.5, random_effects = "x1", id = "id")
  expect_length(yhat, nrow(split$test))
  expect_true(all(yhat %in% c(0, 1)))
})

test_that("fit_gmerf_small is reproducible given the same seed", {
  # NOTE: ranger has stochasticity; we test reproducibility by setting both set.seed()
  # and passing seed= to the fitting function. Keep iterations bounded for runtime.
  df <- sim_data_gmert(G = 6, n_i = 20, seed = 7)

  set.seed(999)
  fit1 <- fit_gmerf_small(
    df = df,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 200,
    max_iter_out = 15,
    tol = 1e-4,
    ntrees = 50,
    seed = 999,
    num.threads = 1
  )

  set.seed(999)
  fit2 <- fit_gmerf_small(
    df = df,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 200,
    max_iter_out = 15,
    tol = 1e-4,
    ntrees = 50,
    seed = 999,
    num.threads = 1
  )

  # core fitted quantities should match exactly (or extremely closely) under same seed
  expect_equal(as.numeric(fit1$sigma2), as.numeric(fit2$sigma2), tolerance = 1e-12)
  expect_equal(fit1$D, fit2$D, tolerance = 1e-12)
  expect_equal(fit1$b, fit2$b, tolerance = 1e-12)

  # predictions should match exactly too
  yhat1 <- predict_gmerf(fit1, new_df = df, thr = 0.5, random_effects = "x1", id = "id")
  yhat2 <- predict_gmerf(fit2, new_df = df, thr = 0.5, random_effects = "x1", id = "id")
  expect_identical(yhat1, yhat2)
})

test_that("fit_gmerf and fit_gmerf_small produce broadly similar fitted values", {
  # NOTE: These models are computationally expensive (outer/inner loops + random forest),
  # so we cap outer iterations to <= 20 to keep test runtime reasonable.
  # The "small" variant uses A^{-1}-based updates and follows a different numerical path,
  # so we do NOT expect exact equality—only broadly similar results.
  #
  # Because the fixed-effect learner is stochastic (ranger), we:
  #  - set.seed()
  #  - pass the same `seed` into both fits
  #  - increase inner iterations to stabilize EM-like updates.

  df <- sim_data_gmert(G = 10, n_i = 30, seed = 101)
  split <- split_gmert_data(df, train_prop = 0.7, seed = 101)

  set.seed(202)
  fit_std <- fit_gmerf(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 500,
    max_iter_out = 20,
    tol = 1e-4,
    ntrees = 100,         # keep moderate for tests
    min_node_size = 5,
    mtry = NULL,
    max.depth = NULL,
    seed = 202
  )

  set.seed(202)
  fit_small <- fit_gmerf_small(
    df = split$train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 500,
    max_iter_out = 20,
    tol = 1e-4,
    ntrees = 100,         # same RF settings as standard
    min_node_size = 5,
    mtry = NULL,
    max.depth = NULL,
    seed = 202,
    num.threads = 1
  )

  # Compare a few core quantities on a *relative* basis
  # (1) residual variance should be close
  expect_equal(as.numeric(fit_std$sigma2), as.numeric(fit_small$sigma2), tolerance = 0.002)

  # (2) random-effects covariance D should be close
  D_std <- fit_std$D
  D_sml <- fit_small$D
  rel_D <- norm(D_std - D_sml, type = "F") / (norm(D_std, type = "F") + 1e-12)
  expect_lt(rel_D, 0.0001)

  # (3) predictions on the test set should be similar (high agreement rate)
  yhat_std <- predict_gmerf(fit_std, new_df = split$test, thr = 0.5, random_effects = "x1", id = "id")
  yhat_sml <- predict_gmerf(fit_small, new_df = split$test, thr = 0.5, random_effects = "x1", id = "id")

  agree <- mean(yhat_std == yhat_sml)
  expect_gt(agree, 0.99)

  # (4) optional: predicted probabilities implied by eta should correlate (if you ever expose them)
  # Here we only have class predictions, so agreement is the main check.
})


test_that("fit_gmerf_small is faster than fit_gmerf", {
  skip_on_cran()

  # keep models small so timings are reasonable but measurable
  df <- sim_data_gmert(G = 8, n_i = 15, seed = 1234)

  set.seed(42)
  t_std <- system.time({
    fit_std_save <- fit_gmerf(
      df = df,
      id = "id",
      target = "y",
      random_effects = "x1",
      max_iter_inn = 100,
      max_iter_out = 100,
      tol = 1e-4,
      ntrees = 20,
      min_node_size = 5,
      mtry = NULL,
      max.depth = NULL,
      seed = 42
        )
  })["elapsed"]

  set.seed(42)
  t_small <- system.time({
    fit_std_nosave <- fit_gmerf_small(
      df = df,
      id = "id",
      target = "y",
      random_effects = "x1",
      max_iter_inn = 100,
      max_iter_out = 100,
      tol = 1e-4,
      ntrees = 20,
      min_node_size = 5,
      mtry = NULL,
      max.depth = NULL,
      seed = 42,
      num.threads = NULL
    )
  })["elapsed"]



  expect_true(is.numeric(t_small) && is.numeric(t_std))

  # small variant should be faster
  expect_lt(as.numeric(t_small), as.numeric(t_std))

})

