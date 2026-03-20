# Predict binary outcomes from a fitted GMERT model

Generates cluster-aware predictions for new data using a fitted
generalized mixed effects regression tree (GMERT) model. Predictions
combine the fixed-effects regression tree component with estimated
cluster-specific random effects when available.

## Usage

``` r
predict_gmert(fit, new_df, random_effect = "x1", thr = 0.5, id = "id")
```

## Arguments

- fit:

  A fitted GMERT model object returned by [`fit_gmert()`](fit_gmert.md).

- new_df:

  A data.frame containing the predictors and cluster identifiers.

- random_effect:

  Character. Name of the random-effects covariate (default: "x1").

- thr:

  Numeric. Classification threshold for predicted probabilities.

- id:

  Character. Name of the cluster identifier column (default: "id").

## Value

Numeric vector of binary predictions (0 or 1).

## References

Hajjem, A., Larocque, D., & Bellavance, F. (2017). *Generalized mixed
effects regression trees*. Statistics & Probability Letters, 126,
114–118.
[doi:10.1016/j.spl.2017.02.033](https://doi.org/10.1016/j.spl.2017.02.033)

## See also

[`fit_gmert()`](fit_gmert.md) for the genral GMERT fitting function.

Other gmert: [`fit_gmerf_small()`](fit_gmerf_small.md),
[`fit_gmert()`](fit_gmert.md),
[`fit_gmert_small()`](fit_gmert_small.md),
[`sim_data_gmert()`](sim_data_gmert.md),
[`split_gmert_data()`](split_gmert_data.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- sim_data_gmert()
split <- split_gmert_data(data)
fit <- fit_gmert(split$train, , max_iter_out = 50)
preds <- predict_gmert(fit, split$test)
table(preds, split$test$y)
} # }
```
