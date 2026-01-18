# Update residual variance using precomputed \\A_i^{-1}\\

Updates the residual variance parameter \\\sigma^2\\ using cluster-wise
quantities and the precomputed matrices \\A_i^{-1}\\ (typically obtained
from [`Ajnv_fun`](Ajnv_fun.md)). This is the “small” (optimized)
counterpart of [`sigma_fun`](sigma_fun.md), avoiding direct
formation/inversion of \\V_i\\.

## Usage

``` r
sigma_fun_small(N, G, idx, Z, w, y_t, fhat, b, Ainv, sigma2)
```

## Arguments

- N:

  Integer. Total number of observations.

- G:

  Integer. Number of clusters.

- idx:

  List of length G. Each element contains the row indices for a given
  cluster.

- Z:

  Numeric matrix (N x q). Random-effects design matrix.

- w:

  Numeric vector of length N. Observation-level weights.

- y_t:

  Numeric vector of length N. Pseudo-response (working response).

- fhat:

  Numeric vector of length N. Current fitted fixed-part predictions.

- b:

  Numeric matrix (G x q). Current random effects (one row per cluster).

- Ainv:

  List of length G. Each element is a (q x q) matrix representing
  \\A_i^{-1}\\ for the corresponding cluster.

- sigma2:

  Numeric scalar. Current residual variance value used in the update.

## Value

A numeric scalar giving the updated residual variance estimate.

## Details

The update aggregates a weighted residual sum of squares across clusters
and includes a trace adjustment term involving
\\\mathrm{tr}(V_i^{-1})\\, computed efficiently via \\A_i^{-1}\\.
