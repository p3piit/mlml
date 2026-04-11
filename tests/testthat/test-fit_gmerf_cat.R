test_that("fit_gmerf_cat returns the expected categorical object structure", {
  df <- gen_sim_cat(
    n = 180,
    n_groups = 12,
    n_vars = 6,
    n_random = 1,
    K = 3,
    seed = 1
  )

  fit <- fit_gmerf_cat(
    df = df,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 8,
    max_iter_out = 8,
    tol = 1e-4,
    ntrees = 50,
    min_node_size = 5,
    sanity_checks = FALSE
  )

  expect_type(fit, "list")
  expect_true(all(c(
    "forests", "b", "D", "sigma2", "mu", "K", "classes",
    "converged_in", "converged_out", "n_iter",
    "train_ids", "tol", "d_eta"
  ) %in% names(fit)))

  expect_true(is.list(fit$forests))
  expect_equal(length(fit$forests), fit$K - 1)
  expect_true(all(vapply(fit$forests, inherits, logical(1), what = "ranger")))

  expect_true(is.matrix(fit$b))
  expect_equal(nrow(fit$b), length(unique(df$id)))

  expect_true(is.matrix(fit$D))
  expect_true(is.numeric(fit$sigma2))
  expect_length(fit$sigma2, 1)
  expect_true(is.finite(fit$sigma2))

  expect_true(is.matrix(fit$mu))
  expect_equal(nrow(fit$mu), nrow(df))
  expect_equal(ncol(fit$mu), fit$K)
  expect_true(all(is.finite(fit$mu)))
  expect_true(all(fit$mu >= 0 & fit$mu <= 1))
  expect_true(all(abs(rowSums(fit$mu) - 1) < 1e-6))

  expect_true(is.numeric(fit$K))
  expect_length(fit$K, 1)
  expect_true(is.character(fit$classes))
  expect_equal(length(fit$classes), fit$K)

  expect_true(is.logical(fit$converged_out))
  expect_length(fit$converged_out, 1)
  expect_true(is.logical(fit$converged_in))

  expect_true(is.numeric(fit$tol))
  expect_length(fit$tol, 1)

  expect_equal(length(fit$train_ids), length(unique(df$id)))
  expect_true(is.numeric(fit$d_eta))
})

test_that("fit_gmerf_cat returns correct stacked dimensions", {
  df <- gen_sim_cat(
    n = 150,
    n_groups = 10,
    n_vars = 5,
    n_random = 1,
    K = 4,
    seed = 2
  )

  fit <- fit_gmerf_cat(
    df = df,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 6,
    max_iter_out = 6,
    tol = 1e-4,
    ntrees = 50
  )

  K1 <- fit$K - 1
  q <- 2  # intercept + x1

  expect_equal(ncol(fit$b), K1 * q)
  expect_equal(dim(fit$D), c(K1 * q, K1 * q))
  expect_equal(length(fit$forests), K1)
})

test_that("fit_gmerf_cat runs end-to-end and predicts classes without error", {
  df <- gen_sim_cat(
    n = 200,
    n_groups = 10,
    n_vars = 6,
    n_random = 1,
    K = 3,
    seed = 3
  )

  fit <- fit_gmerf_cat(
    df = df,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 8,
    max_iter_out = 8,
    tol = 1e-4,
    ntrees = 50,
    sanity_checks = FALSE
  )

  pred <- predict_gmerf_cat(
    fit = fit,
    new_df = df,
    random_effect = "x1",
    id = "id"
  )

  expect_equal(length(pred), nrow(df))
  expect_true(all(pred %in% fit$classes))
})

test_that("fit_gmerf_cat is robust to different seeds/data (smoke test)", {
  df <- gen_sim_cat(
    n = 120,
    n_groups = 6,
    n_vars = 5,
    n_random = 1,
    K = 3,
    seed = 99
  )

  expect_error(
    fit_gmerf_cat(
      df = df,
      id = "id",
      target = "y",
      random_effects = "x1",
      max_iter_inn = 5,
      max_iter_out = 5,
      tol = 1e-4,
      ntrees = 30
    ),
    NA
  )
})

test_that("fit_gmerf_cat fits intercept-only random-effects model", {
  df <- gen_sim_cat(
    n = 120,
    n_groups = 6,
    n_vars = 5,
    n_random = 0,
    K = 3,
    seed = 6
  )

  fit <- fit_gmerf_cat(
    df = df,
    id = "id",
    target = "y",
    random_effects = NULL,
    max_iter_inn = 5,
    max_iter_out = 5,
    tol = 1e-4,
    ntrees = 30
  )

  K1 <- fit$K - 1
  q <- 1  # intercept only

  expect_true(is.matrix(fit$b))
  expect_equal(ncol(fit$b), K1 * q)
  expect_equal(dim(fit$D), c(K1 * q, K1 * q))
  expect_true(!any(is.na(fit$b)))
})

test_that("fit_gmerf_cat rejects binary response", {
  df <- gen_sim(
    n = 100,
    n_groups = 5,
    n_vars = 4,
    seed = 10
  )

  expect_error(
    fit_gmerf_cat(
      df = df,
      id = "id",
      target = "y",
      random_effects = "x1",
      max_iter_inn = 2,
      max_iter_out = 2,
      ntrees = 20
    ),
    "K >= 3"
  )
})

test_that("predict_gmerf_cat works with unseen clusters in new data", {
  df_train <- gen_sim_cat(
    n = 150,
    n_groups = 8,
    n_vars = 5,
    n_random = 1,
    K = 3,
    seed = 11
  )

  df_test <- gen_sim_cat(
    n = 60,
    n_groups = 3,
    n_vars = 5,
    n_random = 1,
    K = 3,
    seed = 12
  )

  df_test$id <- factor(as.integer(df_test$id) + 100)

  fit <- fit_gmerf_cat(
    df = df_train,
    id = "id",
    target = "y",
    random_effects = "x1",
    max_iter_inn = 5,
    max_iter_out = 5,
    tol = 1e-4,
    ntrees = 30
  )

  pred <- predict_gmerf_cat(
    fit = fit,
    new_df = df_test,
    random_effect = "x1",
    id = "id"
  )

  expect_equal(length(pred), nrow(df_test))
  expect_true(all(pred %in% fit$classes))
})