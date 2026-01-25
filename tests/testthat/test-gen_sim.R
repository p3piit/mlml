test_that("basic structure and dimensions are correct", {
  df <- gen_sim(n = 200, n_groups = 10, n_vars = 20, n_random = 3, rho = 0.2, seed = 1)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 200)
  expected_cols <- c("id", "y", paste0("x", 1:20))
  expect_true(all(expected_cols %in% names(df)))
  expect_false(anyNA(df[expected_cols]))
  expect_equal(length(levels(df$id)), 10)
})

test_that("results are reproducible with seed", {
  a <- gen_sim(n = 150, n_groups = 5, n_vars = 10, n_random = 2, rho = 0.1, seed = 42)
  b <- gen_sim(n = 150, n_groups = 5, n_vars = 10, n_random = 2, rho = 0.1, seed = 42)
  expect_equal(a, b)
})

test_that("works for high-dimensional predictor sets and predictors are correlated when rho is large", {
  skip_on_cran()
  # high dimensional run (many predictors)
  df_hd <- gen_sim(n = 300, n_groups = 20, n_vars = 200, n_random = 5, rho = 0.8, seed = 2025)
  expect_equal(nrow(df_hd), 300)
  expect_true(all(paste0("x", 1:200) %in% names(df_hd)))
  # check sample correlation on a subset of X columns reflects high rho
  Xsub <- df_hd[, paste0("x", 1:50)]
  cors <- stats::cor(Xsub)
  mean_offdiag_abs <- mean(abs(cors[lower.tri(cors)]))
  expect_gt(mean_offdiag_abs, 0.5) # with rho=0.8 the average off-diagonal should be substantially > 0.5
})

test_that("invalid random_slope_sd length raises error", {
  expect_error(
    gen_sim(n = 100, n_groups = 5, n_vars = 10, n_random = 3,
                       random_slope_sd = c(0.2, 0.3), seed = 10),
    regexp = "random_slope_sd must be scalar or length n_random"
  )
})