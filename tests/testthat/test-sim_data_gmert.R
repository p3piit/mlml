test_that("sim_data_gmert returns a data.frame with expected columns and types", {
  df <- sim_data_gmert(G = 10, n_i = 20, seed = 1)

  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 10 * 20)

  # expected columns
  expect_true(all(c("id", "y", "x1", "x2", "x3") %in% names(df)))

  # types
  expect_true(is.factor(df$id))
  expect_true(is.numeric(df$y) || is.integer(df$y))
  expect_true(is.numeric(df$x1))
  expect_true(is.numeric(df$x2))
  expect_true(is.numeric(df$x3) || is.integer(df$x3))

  # value constraints
  expect_true(all(df$y %in% c(0, 1)))
  expect_true(all(df$x3 %in% c(0, 1)))

  # x1 range from generator: runif(-2, 2)
  expect_gte(min(df$x1), -2)
  expect_lte(max(df$x1),  2)
})

test_that("sim_data_gmert returns balanced cluster sizes", {
  G <- 7
  n_i <- 13
  df <- sim_data_gmert(G = G, n_i = n_i, seed = 2)

  tab <- table(df$id)
  expect_equal(length(tab), G)
  expect_true(all(as.integer(tab) == n_i))
})

test_that("sim_data_gmert is reproducible with the same seed", {
  df1 <- sim_data_gmert(G = 5, n_i = 10, seed = 123)
  df2 <- sim_data_gmert(G = 5, n_i = 10, seed = 123)

  # strict identical is ok here because the function sets seed internally
  expect_identical(df1, df2)
})

test_that("sim_data_gmert changes when the seed changes", {
  df1 <- sim_data_gmert(G = 5, n_i = 10, seed = 123)
  df2 <- sim_data_gmert(G = 5, n_i = 10, seed = 124)

  # Not guaranteed that everything differs, but extremely likely.
  # We check at least one key column differs.
  expect_false(identical(df1$y, df2$y) && identical(df1$x1, df2$x1))
})

test_that("sim_data_gmert handles edge cases for parameters", {
  # Small but valid sizes
  df <- sim_data_gmert(G = 2, n_i = 1, seed = 1)
  expect_equal(nrow(df), 2)

  # Correlation extremes should not crash in covariance build
  df_rho1  <- sim_data_gmert(G = 3, n_i = 5, rho =  0.99, seed = 1)
  df_rhom1 <- sim_data_gmert(G = 3, n_i = 5, rho = -0.99, seed = 1)

  expect_equal(nrow(df_rho1),  15)
  expect_equal(nrow(df_rhom1), 15)
})

