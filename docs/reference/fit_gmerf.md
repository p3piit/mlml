# Fit a Generalized Mixed-Effects Random Forest model (GMERF)

Fits a Generalized Mixed-Effects Random Forest (GMERF) model for
clustered binary outcomes using an iterative PQL (outer loop) + EM-like
(inner loop) scheme. The fixed-effects component is estimated via a
regression random forest, while the random-effects component is updated
through closed-form updates for \\b\\, \\D\\, and \\\sigma^2\\.

## Usage

``` r
fit_gmerf(
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
  seed = 1234
)
```

## Arguments

- df:

  A data.frame containing the clustered data.

- id:

  Character. Name of the cluster identifier column in `df`.

- target:

  Character. Name of the binary response column in `df`.

- random_effects:

  Character. Name of the covariate used for the random slope; the
  random-effects design is `[1, random_effects]` (random intercept +
  random slope).

- max_iter_inn:

  Integer. Maximum number of inner-loop (EM-like) iterations.

- max_iter_out:

  Integer. Maximum number of outer-loop (PQL) iterations.

- tol:

  Numeric. Convergence tolerance used for both loops.

- ntrees:

  Integer. Number of trees for the random forest.

- mtry:

  Integer or NULL. Number of predictors tried at each split in the
  random forest. If `NULL`, ranger defaults are used.

- min_node_size:

  Integer. Minimum terminal node size for the random forest.

- max.depth:

  Integer or NULL. Maximum depth of trees in the random forest; `NULL`
  means unlimited.

- seed:

  Integer. Random seed passed to the random forest fitting routine.

## Value

A list containing:

- `forest`: fitted ranger forest object (`rf$forest`)

- `b`: estimated random effects matrix (G x q)

- `D`: estimated random-effects covariance matrix (q x q)

- `sigma2`: estimated residual variance

- `mu`: fitted conditional means (length N)

- `converged_in`: logical vector of inner-loop convergence flags (one
  per outer iteration)

- `converged_out`: logical flag for outer-loop convergence

- `n_iter`: number of inner-loop iterations performed in the last outer
  step

- `train_ids`: cluster identifiers used in training (copied from
  `df[[id]]`)

- `gll`: trace of generalized log-likelihood values (for diagnostics)

- `tol`: convergence tolerance used

## Details

The implementation follows the generalized mixed-effects regression tree
/ forest framework, with the fixed-effects learner replaced by a random
forest (see References).

This routine relies on helper functions [`Vi_fun()`](Vi_fun.md),
[`b_fun()`](b_fun.md), [`sigma_fun()`](sigma_fun.md),
[`D_fun()`](D_fun.md), and [`gll_fun()`](gll_fun.md) for mixed-effects
updates, and uses
[`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md)
for the fixed-effects regression learner.

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

[`predict_gmerf`](predict_gmerf.md) for prediction on new data.

[`fit_gmert`](fit_gmert.md) for the tree-based analogue.

Other gmerf: [`predict_gmerf()`](predict_gmerf.md)

## Examples

``` r
library(ranger)
library(MASS)
# Simulate clustered binary data
df <- sim_data_gmert(G = 20, n_i = 30, seed = 1)

# Fit GMERF
fit <- fit_gmerf(df, ntrees = 200, seed = 1, max_iter_out = 50)
#> WARNING: the PQL algorithm did not converge in 50 iterations.

# Inspect convergence
fit$converged_out
#> [1] FALSE
```
