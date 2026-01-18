# Update random-effects covariance matrix D

Internal helper for GMERT/GMERF. Updates the random-effects covariance
matrix \\D\\ using the closed-form expression based on the empirical
second moments of the random effects and a shrinkage correction term
derived from the conditional covariance.

## Usage

``` r
D_fun(G, idx, b, Z, D, Vi, w)
```

## Arguments

- G:

  Integer. Number of clusters.

- idx:

  List of length `G`. Row indices for each cluster.

- b:

  Numeric matrix (G x q). Current random-effects estimates.

- Z:

  Numeric matrix (N x q). Random-effects design matrix.

- D:

  Numeric matrix (q x q). Current random-effects covariance matrix.

- Vi:

  List of length `G`. Cluster-specific covariance matrices \\V_i\\.

- w:

  Numeric vector length N. Working weights.

## Value

Numeric matrix (q x q). Updated estimate of the random-effects
covariance matrix \\D\\.
