# Train/test split by cluster for GMERT data

Splits a clustered dataset into training and test sets while preserving
the cluster structure. Each cluster contributes the same proportion of
observations to the training and test sets.

## Usage

``` r
split_gmert_data(df, train_prop = 0.7, seed = 123)
```

## Arguments

- df:

  A data.frame containing clustered observations. Must include a column
  named `id`.

- train_prop:

  Numeric. Proportion of observations per cluster to use for training.

- seed:

  Integer. Random seed for reproducibility.

## Value

A list with two elements:

- train:

  Training data.frame

- test:

  Test data.frame

## See also

Other gmert: [`fit_gmerf_small()`](fit_gmerf_small.md),
[`fit_gmert()`](fit_gmert.md),
[`fit_gmert_small()`](fit_gmert_small.md),
[`predict_gmert()`](predict_gmert.md),
[`sim_data_gmert()`](sim_data_gmert.md)

## Examples

``` r
data <- sim_data_gmert()
split <- split_gmert_data(data)
nrow(split$train)
#> [1] 1400
nrow(split$test)
#> [1] 600
```
