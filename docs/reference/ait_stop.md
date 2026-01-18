# Aitken acceleration convergence check

Internal helper for GMERT/GMERF. Applies Aitken acceleration to the
generalized log-likelihood sequence in order to estimate the asymptotic
value and assess convergence based on a relative error criterion.

## Usage

``` r
ait_stop(gll, tol, it)
```

## Arguments

- gll:

  Numeric vector. History of generalized log-likelihood values (must
  contain at least three consecutive iterations).

- tol:

  Numeric scalar. Tolerance for the relative error to the extrapolated
  asymptote.

- it:

  Integer. Current iteration index (referring to the position in `gll`).

## Value

Logical. `TRUE` if convergence is detected, `FALSE` otherwise.
