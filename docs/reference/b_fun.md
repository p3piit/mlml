# Update cluster-level random effects

Internal update step for GMERT/GMERF. Computes the cluster-specific
random effects \\b_g\\ given current covariance components and working
quantities.

## Usage

``` r
b_fun(G, Z, Vi, D, idx, y_t, fhat, w)
```

## Arguments

- G:

  Integer. Number of clusters.

- Z:

  Numeric matrix (N x q). Random-effects design matrix (stacked by
  cluster).

- Vi:

  List of length `G`. Each element is an \\n_g \times n_g\\ covariance
  matrix for cluster `g`.

- D:

  Numeric matrix (q x q). Random-effects covariance matrix.

- idx:

  List of length `G`. Row indices for each cluster in `Z` / `y_t` /
  `fhat` / `w`.

- y_t:

  Numeric vector length N. Working (pseudo-)response.

- fhat:

  Numeric vector length N. Fixed-part prediction.

- w:

  Numeric vector length N. Working weights (typically \\\mu (1-\mu)\\).

## Value

A numeric matrix (G x q) of updated random effects.
