# Build cluster-specific covariance matrices V_i

Internal helper for GMERT/GMERF. For each cluster \\i\\, constructs
\\V_i = (W_i^{1/2} Z_i) D (W_i^{1/2} Z_i)^T + \sigma^2 I\_{n_i}\\, where
\\W_i\\ is the diagonal working-weight matrix.

## Usage

``` r
Vi_fun(Z, D, s2, G, idx, w)
```

## Arguments

- Z:

  Numeric matrix (N x q). Random-effects design matrix (stacked by
  cluster).

- D:

  Numeric matrix (q x q). Random-effects covariance matrix.

- s2:

  Numeric scalar. Residual variance on the working scale.

- G:

  Integer. Number of clusters.

- idx:

  List of length `G`. Row indices for each cluster.

- w:

  Numeric vector length N. Working weights.

## Value

A list of length `G`. Each element is an \\n_i \times n_i\\ covariance
matrix \\V_i\\.
