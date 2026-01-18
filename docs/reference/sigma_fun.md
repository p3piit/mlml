# Update residual variance \\\sigma^2\\

Internal helper for GMERT/GMERF. Updates the residual variance parameter
\\\sigma^2\\ using the closed-form expression derived from the
generalized log-likelihood, aggregating cluster-wise weighted residuals.

## Usage

``` r
sigma_fun(N, G, idx, b, y, Z, D, Vi, s2, fhat, w)
```

## Arguments

- N:

  Integer. Total number of observations.

- G:

  Integer. Number of clusters.

- idx:

  List of length `G`. Row indices for each cluster.

- b:

  Numeric matrix (G x q). Current random-effects estimates.

- y:

  Numeric vector length N. Response or pseudo-response.

- Z:

  Numeric matrix (N x q). Random-effects design matrix.

- D:

  Numeric matrix (q x q). Random-effects covariance matrix.

- Vi:

  List of length `G`. Cluster-specific covariance matrices \\V_i\\.

- s2:

  Numeric scalar. Current residual variance.

- fhat:

  Numeric vector length N. Fixed-effects predictions.

- w:

  Numeric vector length N. Working weights.

## Value

Numeric scalar. Updated estimate of \\\sigma^2\\.
