# ait_stop: Aitken acceleration based stopping rule
#'
#' Aitken acceleration convergence check
#'
#' Internal helper for GMERT/GMERF. Applies Aitken acceleration to the
#' generalized log-likelihood sequence in order to estimate the asymptotic
#' value and assess convergence based on a relative error criterion.
#'
#' @param gll Numeric vector. History of generalized log-likelihood values
#'   (must contain at least three consecutive iterations).
#' @param tol Numeric scalar. Tolerance for the relative error to the
#'   extrapolated asymptote.
#' @param it Integer. Current iteration index (referring to the position in `gll`).
#'
#' @return Logical. `TRUE` if convergence is detected, `FALSE` otherwise.
#'
#' @keywords internal
ait_stop <- function(gll,  #gll : numeric vector with GLL history (needs last three values)
                     tol,  #tol : tolerance for relative error to asymptote
                     it    #it  : current iteration index
){
  converged <- FALSE
  l1 <- gll[it + 2] - gll[it + 1]                                      # L_t
  l2 <- gll[it + 1] - gll[it]                                          # L_{t-1}
  a <- l1 / l2                                                         # Aitken ratio
  if (is.finite(a) && a > 0 && a < 1) {                                # only meaningful if 0<a<1
    L_inf <- gll[it + 1] + l1 / (1 - a)                                # extrapolate asymptote
    err <- abs(L_inf - gll[it + 2]) / (abs(L_inf) + 1e-12)             # relative error to asymptote
    if (err < tol) converged <- TRUE                                   # declare convergence
  }
  converged                                                            # return flag
}
