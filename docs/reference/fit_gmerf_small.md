# Fit a generalized mixed-effects random forest (GMERF) using A-inverse updates

Fits a generalized mixed-effects model where the fixed-effects component
is learned via a random forest and the random effects account for
clustered/longitudinal dependence. Estimation proceeds via a PQL outer
loop and an EM-like inner loop, using precomputed cluster-specific
\\A_i^{-1}\\ matrices to update random effects and variance components
efficiently.

## Usage

``` r
fit_gmerf_small(
  df,
  id = "id",
  target = "y",
  random_effects = "x1",
  max_iter_inn = 1000,
  max_iter_out = 1000,
  tol = 1e-06,
  ntrees = 500,
  mtry = NULL,
  min_node_size = 5,
  max.depth = NULL,
  seed = 1234,
  num.threads = NULL
)
```

## Arguments

- df:

  A data.frame containing the response, cluster identifier,
  random-effect covariate(s), and fixed-effect predictors.

- id:

  Character scalar. Name of the clustering/grouping variable in `df`.

- target:

  Character scalar. Name of the binary response variable in `df`.

- random_effects:

  Character scalar. Name(s) of the variable(s) used for the
  random-effects design (intercept is added internally).

- max_iter_inn:

  Integer. Maximum number of EM-like iterations (inner loop).

- max_iter_out:

  Integer. Maximum number of PQL iterations (outer loop).

- tol:

  Numeric. Convergence tolerance used for both inner and outer loops.

- ntrees:

  Integer. Number of trees for the random forest.

- mtry:

  Integer or `NULL`. Number of candidate features tried at each split
  (passed to
  [`ranger::ranger`](http://imbs-hl.github.io/ranger/reference/ranger.md);
  `NULL` uses `ranger`'s default).

- min_node_size:

  Integer. Minimum terminal node size for the forest.

- max.depth:

  Integer or `NULL`. Maximum tree depth for the forest (`NULL` =
  unlimited).

- seed:

  Integer. Random seed for reproducibility.

- num.threads:

  Integer or `NULL`. Number of threads used by `ranger` (`NULL` uses all
  available threads).

## Value

A list with components:

- forest:

  Fitted random forest object (stored as `rf$forest`).

- b:

  Matrix of estimated random effects (clusters x q).

- D:

  Estimated random-effects covariance matrix.

- sigma2:

  Estimated residual variance.

- mu:

  Fitted conditional means.

- converged_in:

  Logical vector of inner-loop convergence flags (per outer iteration).

- converged_out:

  Logical flag for outer-loop convergence.

- n_iter:

  Number of inner-loop iterations performed in the last outer iteration.

- train_ids:

  Cluster identifiers used for training.

- gll:

  Generalized log-likelihood trace used for diagnostics.

- tol:

  Tolerance used for convergence checks.

## References

Hajjem, Ahlem, François Bellavance, and Denis Larocque. 2011. “Mixed
Effects Regression Trees for Clustered Data.” *Statistics & Probability
Letters* 81 (4): 451–59.
[doi:10.1016/j.spl.2010.12.003](https://doi.org/10.1016/j.spl.2010.12.003)
.

Hajjem, Ahlem, Denis Larocque, and François Bellavance. 2017.
“Generalized Mixed Effects Regression Trees.” *Statistics & Probability
Letters* 126 (July): 114–18.
[doi:10.1016/j.spl.2017.02.033](https://doi.org/10.1016/j.spl.2017.02.033)
.

Pellagatti, Massimo, Chiara Masci, Francesca Ieva, and Anna M. Paganoni.
2021. “Generalized Mixed-Effects Random Forest: A Flexible Approach to
Predict University Student Dropout.” *Statistical Analysis and Data
Mining: The ASA Data Science Journal* 14 (3): 241–57.
[doi:10.1002/sam.11505](https://doi.org/10.1002/sam.11505) .

## See also

[`fit_gmerf`](fit_gmerf.md) for the baseline GMERF fitting routine;
[`fit_gmert_small`](fit_gmert_small.md) for the tree-based \\A^{-1}\\
optimized variant.

Other gmert: [`fit_gmert()`](fit_gmert.md),
[`fit_gmert_small()`](fit_gmert_small.md),
[`predict_gmert()`](predict_gmert.md),
[`sim_data_gmert()`](sim_data_gmert.md),
[`split_gmert_data()`](split_gmert_data.md)

## Examples

``` r
# Simulate clustered binary data
df <- sim_data_gmert(G = 10, n_i = 30, seed = 1)

# Fit the optimized GMERF (A^{-1}-based updates)
fit <- fit_gmerf_small(df, id = "id", target = "y", random_effects = "x1",
                       max_iter_inn = 50, max_iter_out = 50, tol = 1e-4,
                       ntrees = 100, seed = 1)
#> WARNING: the PQL algorithm did not converge in 50 iterations.

# Predicted classes on training data
pred <- predict_gmerf(fit, new_df = df, random_effects = "x1", id = "id")
```
