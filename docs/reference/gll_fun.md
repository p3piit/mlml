# Compute generalized log-likelihood (up to an additive constant)

Internal helper for GMERT/GMERF. Computes the generalized log-likelihood
criterion used for monitoring and convergence checks. The returned value
is defined up to an additive constant (i.e., it is suitable for
comparing iterations but not as an absolute likelihood value).

## Usage

``` r
gll_fun(idx, b, y, Z, D, s2, fhat, w)
```

## Arguments

- idx:

  List of row indices per cluster.

- b:

  Numeric matrix (G x q). Random-effects estimates.

- y:

  Numeric vector length N. Response or pseudo-response.

- Z:

  Numeric matrix (N x q). Random-effects design matrix.

- D:

  Numeric matrix (q x q). Random-effects covariance matrix.

- s2:

  Numeric scalar. Residual variance.

- fhat:

  Numeric vector length N. Fixed-effects predictions.

- w:

  Numeric vector length N. Working weights.

## Value

Numeric scalar. Generalized log-likelihood value (up to an additive
constant).
