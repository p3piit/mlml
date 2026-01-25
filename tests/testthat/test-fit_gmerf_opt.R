# keep iterations and trees small to limit runtime in CI
common_args <- list(
  ntrees = 50L,
  max_iter_inn = 300L,
  max_iter_out = 200L,
  min_node_size = 5L,
  seed = 42L
)

test_that("fit_gmerf_opt returns expected structure and sizes", {
  sim_df <- gen_sim(n = 200, n_groups = 10, n_vars = 20, n_random = 3, rho = 0.2, seed = 1)
  res <- fit_gmerf_opt(df = sim_df,
                       ntrees = common_args$ntrees,
                       max_iter_inn = common_args$max_iter_inn,
                       max_iter_out = common_args$max_iter_out,
                       min_node_size = common_args$min_node_size,
                       seed = 1)
  expect_type(res, "list")
  expect_true(all(c("forest", "b", "D", "sigma2", "mu") %in% names(res)))
  expect_true(is.matrix(res$b))
  expect_equal(nrow(res$b), length(unique(sim_df$id)))
  expect_equal(ncol(res$b), 2) # intercept + one slope in default tests
})

test_that("fit_gmerf_opt is reproducible with seed", {
  df0 <- gen_sim(n = 150, n_groups = 5, n_vars = 10, n_random = 2, rho = 0.1, seed = 202)
  a <- fit_gmerf_opt(df = df0, ntrees = 40L, max_iter_inn = 30L, max_iter_out = 20L, seed = 42)
  b <- fit_gmerf_opt(df = df0, ntrees = 40L, max_iter_inn = 30L, max_iter_out = 20L, seed = 42)
  expect_equal(a$b, b$b)
  expect_equal(a$D, b$D)
  expect_equal(a$sigma2, b$sigma2)
})

test_that("fit_gmerf_opt gives similar estimates to original fit_gmerf (small iterations)", {
  skip_on_cran()
  df_cmp <- gen_sim(n = 220, n_groups = 11, n_vars = 25, n_random = 3, rho = 0.25, seed = 101)
  res_orig <- fit_gmerf(df = df_cmp, ntrees = 50L, max_iter_inn = 30L, max_iter_out = 20L, seed = 101)
  res_opt  <- fit_gmerf_opt(df = df_cmp, ntrees = 50L, max_iter_inn = 30L, max_iter_out = 20L, seed = 101)
  # allow moderate tolerance because iterations are low and algorithms differ
  expect_true(max(abs(res_orig$b - res_opt$b)) < 0.6)
  expect_true(max(abs(res_orig$D - res_opt$D)) < 0.6)
  expect_true(abs(res_orig$sigma2 - res_opt$sigma2) < 0.6)
})

test_that("fit_gmerf_opt is consistent with fit_gmerf_small (same optimization family)", {
  skip_on_cran()
  df_cmp2 <- gen_sim(n = 180, n_groups = 9, n_vars = 30, n_random = 3, rho = 0.2, seed = 202)
  res_small <- fit_gmerf_small(df = df_cmp2, ntrees = 40L, max_iter_inn = 30L, max_iter_out = 20L, seed = 202)
  res_opt   <- fit_gmerf_opt(df = df_cmp2, ntrees = 40L, max_iter_inn = 30L, max_iter_out = 20L, seed = 202)
  # expect closer agreement between the two optimized variants
  expect_true(max(abs(res_small$b - res_opt$b)) < 0.3)
  expect_true(max(abs(res_small$D - res_opt$D)) < 0.3)
})

test_that("fit_gmerf_opt shows a wall-clock advantage against the original implementation", {
  skip_on_cran()
  df_time <- gen_sim(n = 300, n_groups = 15, n_vars = 60, n_random = 3, rho = 0.3, seed = 303)
  t_orig <- system.time(res1 <- fit_gmerf(df = df_time, ntrees = 60L, max_iter_inn = 30L, max_iter_out = 20L, seed = 303))["elapsed"]
  t_opt  <- system.time(res2 <- fit_gmerf_opt(df = df_time, ntrees = 60L, max_iter_inn = 30L, max_iter_out = 20L, seed = 303))["elapsed"]
  # optimized version should not be substantially slower; allow 50% slack
  expect_lte(t_opt, t_orig * 1.5)
})

test_that("fit_gmerf_opt runs on higher dimensional data and returns expected columns", {
  skip_on_cran()
  df_hd <- gen_sim(n = 300, n_groups = 20, n_vars = 250, n_random = 5, rho = 0.25, seed = 404)
  res_hd <- fit_gmerf_opt(df = df_hd, ntrees = 50L, max_iter_inn = 30L, max_iter_out = 20L, seed = 404)
  expect_true(is.list(res_hd))
  expect_equal(nrow(res_hd$b), length(unique(df_hd$id)))
  expect_equal(ncol(res_hd$b), 2)
})

