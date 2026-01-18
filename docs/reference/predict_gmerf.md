# Predict class labels from a fitted GMERF model

Generates binary predictions for new clustered data using a fitted
Generalized Mixed-Effects Random Forest (GMERF) model. Predictions
combine the fixed-effects random forest component with cluster-specific
random effects estimated during training.

## Usage

``` r
predict_gmerf(fit, new_df, thr = 0.5, random_effects = "x1", id = "id")
```

## Arguments

- fit:

  A fitted GMERF model object returned by [`fit_gmerf`](fit_gmerf.md).

- new_df:

  A data.frame containing new observations for prediction.

- thr:

  Numeric. Classification threshold for converting predicted
  probabilities into class labels (default: 0.5).

- random_effects:

  Character. Name of the covariate used for the random slope in the
  random-effects design.

- id:

  Character. Name of the cluster identifier column in `new_df`.

## Value

A numeric vector of predicted class labels (0/1).

## Details

For clusters observed in the training data, the corresponding random
effects are added to the fixed-effects prediction. For unseen clusters,
only the fixed-effects component is used.

The linear predictor is computed as: \$\$\eta = \hat f(X) + Z_i \hat
b_i\$\$ where \\\hat f(X)\\ is the random forest prediction and \\Z_i
\hat b_i\\ is the random-effects contribution for cluster \\i\\, when
available.

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

[`fit_gmerf`](fit_gmerf.md) for model fitting.

Other gmerf: [`fit_gmerf()`](fit_gmerf.md)

## Examples

``` r
# Simulate clustered data
df <- sim_data_gmert(G = 20, n_i = 30, seed = 1)

# Fit GMERF
fit <- fit_gmerf(df, ntrees = 200, seed = 1, , max_iter_out = 50)
#> WARNING: the PQL algorithm did not converge in 50 iterations.

# Predict on training data
yhat <- predict_gmerf(fit, df)
table(yhat, df$y)
#>     
#> yhat   0   1
#>    0 191  37
#>    1  52 320
```
