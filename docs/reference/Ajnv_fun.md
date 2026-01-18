# Precompute cluster-specific inverse matrices for random-effects updates

Computes the inverse matrices \\A_j^{-1}\\ for each cluster \\j\\, where
\\A_j = D^{-1} + (1 / \sigma^2) Z_j^T W_j Z_j\\. These matrices are used
to speed up the closed-form updates of random effects and variance
components in the small-sample GMERT / GMERF algorithms.

## Usage

``` r
Ajnv_fun(G, Z, w, D, sigma2, idx)
```

## Arguments

- G:

  Integer. Number of clusters.

- Z:

  Numeric matrix (N x q). Random-effects design matrix.

- w:

  Numeric vector of length N. Observation-level weights.

- D:

  Numeric matrix (q x q). Current random-effects covariance matrix.

- sigma2:

  Numeric scalar. Current residual variance.

- idx:

  List of length G. Each element contains the row indices of `Z`
  corresponding to a given cluster.

## Value

A list of length G, where each element is a (q x q) matrix containing
\\A_j^{-1}\\ for the corresponding cluster.

## Details

The matrices returned by this function are used in alternative update
formulas for the random effects \\b_j\\, the residual variance
\\\sigma^2\\, and the random-effects covariance matrix \\D\\, reducing
computational cost by avoiding repeated matrix inversions.
