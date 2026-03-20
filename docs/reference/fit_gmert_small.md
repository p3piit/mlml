# Fit a generalized mixed-effects regression tree (GMERT) using A-inverse updates

Fits a generalized mixed-effects regression tree model for clustered
binary data, using an optimized set of updates based on precomputing
\\A_i^{-1}\\ to speed up estimation of the random effects covariance and
variance components.

## Usage

``` r
fit_gmert_small(
  df,
  id = "id",
  target = "y",
  random_effects = "x1",
  max_iter_inn = 1000,
  max_iter_out = 1000,
  tol = 1e-06,
  cp = 0,
  minsplit = 50,
  minbucket = 20,
  maxdepth = 5,
  xval = 10
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

- cp:

  Numeric. Complexity parameter passed to `rpart.control()`.

- minsplit:

  Integer. Minimum number of observations required to attempt a split.

- minbucket:

  Integer. Minimum number of observations in any terminal node.

- maxdepth:

  Integer. Maximum depth of the fitted tree.

- xval:

  Integer. Number of cross-validation folds in `rpart`.

## Value

A list with components:

- tree:

  Fitted `rpart` regression tree for the fixed-effects component.

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

- time:

  Elapsed runtime in seconds.

## Details

This routine follows the penalized quasi-likelihood (PQL) outer loop and
an EM-like inner loop described for GMERT, while replacing direct
\\V_i^{-1}\\ computations with cluster-wise precomputed \\A_i^{-1}\\
matrices.

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

## See also

[`fit_gmert`](fit_gmert.md) for the baseline GMERT fitting routine using
\\V_i\\.

Other gmert: [`fit_gmerf_small()`](fit_gmerf_small.md),
[`fit_gmert()`](fit_gmert.md), [`predict_gmert()`](predict_gmert.md),
[`sim_data_gmert()`](sim_data_gmert.md),
[`split_gmert_data()`](split_gmert_data.md)

## Examples

``` r
# Simulate clustered binary data
df <- sim_data_gmert(G = 10, n_i = 30, seed = 1)

# Fit the optimized GMERT (A^{-1}-based updates)
fit <- fit_gmert_small(df, id = "id", target = "y", random_effects = "x1",
                       max_iter_inn = 50, max_iter_out = 50, tol = 1e-4)
#> WARNING: the PQL algorithm did not converge in 50 iterations.

# Predicted classes on training data
pred <- predict_gmert(fit, new_df = df, random_effect = "x1", id = "id")
```
