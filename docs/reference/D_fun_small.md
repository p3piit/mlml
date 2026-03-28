# Update random-effects covariance matrix using precomputed \\A_i^{-1}\\

Updates the random-effects covariance matrix \\D\\ by combining the
outer products of the current random effects \\b_i\\ with the
corresponding precomputed inverse matrices \\A_i^{-1}\\. This function
is the optimized (“small”) counterpart of `D_fun`, avoiding direct
dependence on \\V_i\\.

## Usage

``` r
D_fun_small(G, b, Ainv)
```

## Arguments

- G:

  Integer. Number of clusters.

- b:

  Numeric matrix (G x q). Current random effects, one row per cluster.

- Ainv:

  List of length G. Each element is a (q x q) matrix representing
  \\A_i^{-1}\\ for the corresponding cluster.

## Value

A numeric (q x q) matrix giving the updated estimate of the
random-effects covariance matrix \\D\\.

## Details

The update is computed as \$\$D = \frac{1}{G} \sum\_{i=1}^G \left( b_i
b_i^T + A_i^{-1} \right),\$\$ which combines the empirical covariance of
the random effects with a shrinkage term derived from the inverse
information matrices.
