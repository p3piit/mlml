# Compute the F1-score for the majority or minority class

Computes the F1-score from a 2x2 confusion matrix, either for the
majority class (default) or for the minority class.

## Usage

``` r
f1_fun(table, majority = TRUE)
```

## Arguments

- table:

  A 2x2 confusion matrix (e.g., produced by `table(y, yhat)`), where
  rows correspond to true classes and columns to predicted classes.

- majority:

  Logical. If `TRUE`, compute the F1-score for the majority class; if
  `FALSE`, compute it for the minority class.

## Value

A numeric scalar giving the F1-score for the selected class.

## Examples

``` r
# Example confusion matrix
cm <- matrix(c(50, 10,
               5,  35), nrow = 2, byrow = TRUE)

# F1-score for the majority class
f1_fun(cm, majority = TRUE)
#> [1] 0.8695652

# F1-score for the minority class
f1_fun(cm, majority = FALSE)
#> [1] 0.8235294
```
