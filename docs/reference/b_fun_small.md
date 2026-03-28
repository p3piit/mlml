# Update cluster-specific random effects using precomputed \\A_i^{-1}\\

Computes the updated random effects \\b_i\\ for each cluster \\i\\ using
the precomputed inverse matrices \\A_i^{-1}\\ (typically returned by
[`Ajnv_fun`](Ajnv_fun.md)). This version avoids forming and inverting
\\V_i\\ directly, and is used in the “small” GMERT / GMERF optimization
routines.

## Usage

``` r
b_fun_small(G, Z, w, idx, y_t, fhat, Ainv, sigma2)
```

## Arguments

- G:

  Integer. Number of clusters.

- Z:

  Numeric matrix (N x q). Random-effects design matrix.

- w:

  Numeric vector of length N. Observation-level weights.

- idx:

  List of length G. Each element contains the row indices of `Z` (and
  associated vectors) corresponding to a given cluster.

- y_t:

  Numeric vector of length N. Pseudo-response (working response).

- fhat:

  Numeric vector of length N. Current fitted fixed-part predictions.

- Ainv:

  List of length G. Each element is a (q x q) matrix representing
  \\A_i^{-1}\\ for the corresponding cluster.

- sigma2:

  Numeric scalar. Current residual variance.

## Value

A numeric matrix of dimension (G x q) containing the updated random
effects for each cluster.

## Details

For each cluster \\i\\, this function computes \$\$b_i = A_i^{-1} \cdot
(1/\sigma^2) Z_i^T W_i (y^\*\_i - f_i)\$\$ where \\y^\*\_i\\ is the
pseudo-response restricted to cluster \\i\\, \\f_i\\ is the fixed-part
prediction, and \\W_i\\ is a diagonal matrix of working weights.
