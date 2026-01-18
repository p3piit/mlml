# Fit a generalized mixed-effects regression tree (GMERT)

Fits a generalized mixed-effects regression tree model for clustered
binary outcomes using a PQL outer loop and an EM-like inner loop. The
fixed part is estimated via an `rpart` regression tree on a working
pseudo-response, and the random part is modeled with cluster-specific
random intercepts and slopes (currently implemented as an intercept plus
one random slope).

## Usage

``` r
fit_gmert(
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

  A data frame containing the response, cluster identifier, and
  predictors.

- id:

  Character string. Name of the cluster identifier column in `df`.

- target:

  Character string. Name of the binary response column in `df` (coded
  0/1).

- random_effects:

  Character string. Name of the predictor used for the random slope. The
  random-effects design is built as an intercept plus this column.

- max_iter_inn:

  Integer. Maximum number of EM-like inner-loop iterations.

- max_iter_out:

  Integer. Maximum number of PQL outer-loop iterations.

- tol:

  Numeric scalar. Convergence tolerance used for both inner and outer
  loops.

- cp:

  Numeric scalar. `rpart` complexity parameter (pruning threshold).

- minsplit:

  Integer. Minimum number of observations required to attempt a split in
  `rpart`.

- minbucket:

  Integer. Minimum number of observations in any terminal node in
  `rpart`.

- maxdepth:

  Integer. Maximum depth of the fitted `rpart` tree.

- xval:

  Integer. Number of cross-validation folds used by `rpart`.

## Value

A list with components:

- tree:

  An `rpart` object: the fitted regression tree for the fixed effects.

- b:

  Numeric matrix (G x q) of estimated random effects by cluster.

- D:

  Numeric matrix (q x q): estimated random-effects covariance matrix.

- sigma2:

  Numeric scalar: estimated residual variance on the working scale.

- mu:

  Numeric vector length N: fitted conditional means at convergence.

- converged_in:

  Logical vector: inner-loop convergence flag for each outer iteration.

- converged_out:

  Logical scalar: outer-loop convergence flag.

- n_iter:

  Integer: number of inner iterations performed in the last outer
  iteration.

- train_ids:

  Vector of cluster identifiers from the training data.

- gll:

  Numeric vector: generalized log-likelihood trace used for diagnostics.

- tol:

  Numeric scalar: tolerance used.

## Details

The function returns a list containing the fitted tree, estimated random
effects, covariance parameters, and convergence diagnostics.

This implementation assumes a binary response and uses a logistic link
via a penalized quasi-likelihood (PQL) linearization. The random-effects
structure is currently an intercept plus one slope specified by
`random_effects`.

## References

Hajjem, A., Larocque, D., & Bellavance, F. (2017). *Generalized mixed
effects regression trees*. Statistics & Probability Letters, 126,
114–118.
[doi:10.1016/j.spl.2017.02.033](https://doi.org/10.1016/j.spl.2017.02.033)

## See also

[`predict_gmert()`](predict_gmert.md) for prediction on new data.

Other gmert: [`fit_gmerf_small()`](fit_gmerf_small.md),
[`fit_gmert_small()`](fit_gmert_small.md),
[`predict_gmert()`](predict_gmert.md),
[`sim_data_gmert()`](sim_data_gmert.md),
[`split_gmert_data()`](split_gmert_data.md)

## Examples

``` r
# Simulate clustered binary data
df <- sim_data_gmert(G = 10, n_i = 30, seed = 1)
split <- split_gmert_data(df, train_prop = 0.7, seed = 1)

# Fit GMERT
fit <- fit_gmert(split$train, id = "id", target = "y", random_effects = "x1",
                 max_iter_inn = 50, max_iter_out = 50, tol = 1e-4,
                 maxdepth = 3, minsplit = 20, minbucket = 10)
#> WARNING: the PQL algorithm did not converge in 50 iterations.

# Predict on test set
pred <- predict_gmert(fit, split$test, random_effect = "x1", id = "id")
table(truth = split$test$y, pred = pred)
#>      pred
#> truth  0  1
#>     0 18 20
#>     1 10 42
```
