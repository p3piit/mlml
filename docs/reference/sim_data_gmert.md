# Simulate clustered binary data for GMERT models

Generates synthetic clustered data from a logistic mixed-effects model
with random intercepts and random slopes. The simulated dataset is
suitable for fitting and testing generalized mixed effects regression
tree (GMERT) models.

## Usage

``` r
sim_data_gmert(
  G = 50,
  n_i = 40,
  beta0 = 0.5,
  beta1 = 1.2,
  beta2 = -0.8,
  beta3 = 0.6,
  sigma_b0 = 0.8,
  sigma_b1 = 0.5,
  rho = 0.2,
  seed = 123
)
```

## Arguments

- G:

  Integer. Number of clusters.

- n_i:

  Integer. Number of observations per cluster.

- beta0:

  Numeric. Fixed intercept.

- beta1:

  Numeric. Fixed effect for `x1`.

- beta2:

  Numeric. Fixed effect for `x2`.

- beta3:

  Numeric. Fixed effect for `x3`.

- sigma_b0:

  Numeric. Standard deviation of the random intercept.

- sigma_b1:

  Numeric. Standard deviation of the random slope for `x1`.

- rho:

  Numeric. Correlation between random intercept and random slope.

- seed:

  Integer. Random seed for reproducibility.

## Value

A data.frame with columns `id`, `y`, `x1`, `x2`, and `x3`.

## See also

[`fit_gmert`](fit_gmert.md), [`predict_gmert`](predict_gmert.md)

Other gmert: [`fit_gmerf_small()`](fit_gmerf_small.md),
[`fit_gmert()`](fit_gmert.md),
[`fit_gmert_small()`](fit_gmert_small.md),
[`predict_gmert()`](predict_gmert.md),
[`split_gmert_data()`](split_gmert_data.md)

## Examples

``` r
data <- sim_data_gmert()
head(data)
#>   id y          x1         x2 x3
#> 1  1 1 -1.04509589  0.9159921  1
#> 2  1 1  1.84943575  0.8006224  1
#> 3  1 0  0.40546290 -0.9365690  0
#> 4  1 1  0.06011891 -1.4007874  1
#> 5  1 1 -0.38970663  0.1602775  1
#> 6  1 1  1.52098616 -0.2739624  0
```
