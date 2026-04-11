#'
#' Predict categorical outcomes from a fitted GMERT model
#'
#' Generates cluster-aware predictions for new data using a fitted
#' generalized mixed effects regression tree (GMERT) model. Predictions
#' combine the fixed-effects regression tree component with estimated
#' cluster-specific random effects when available.
#'
#' @param fit A fitted GMERT model object returned by \code{fit_gmert_small()}.
#' @param new_df A data.frame containing the predictors and cluster identifiers.
#' @param random_effect Character. Name of the random-effects covariate
#'   (default: \code{"x1"}).
#' @param id Character. Name of the cluster identifier column
#'   (default: \code{"id"}).
#'
#' @return
#' Returns a vector of predicted classes.
#'
#' @family gmert
#'
#' @export
predict_gmert_cat <- function(fit,
                          new_df,
                          random_effect = "x1",
                          id = "id") {

  
  N_new <- nrow(new_df)
  K <- fit$K
  K1 <- K - 1L

  # fixed part
  fhat <- matrix(0, N_new, K1)
  for (k in seq_len(K1)) {
    fhat[, k] <- as.numeric(predict(fit$trees[[k]], newdata = new_df))
  }

  # random-effects design
  if (is.null(random_effect)) {
  Znew <- matrix(1, nrow = nrow(new_df), ncol = 1)
} else {
  Znew <- as.matrix(cbind(1, new_df[random_effect]))
}
  q <- ncol(Znew)

  # random-effects contribution
  add <- matrix(0, N_new, K1)

  clus_fit <- fit$train_ids
  map <- match(new_df[[id]], clus_fit)
  seen <- !is.na(map)

  if (any(seen)) {
    for (r in which(seen)) {
      g <- map[r]
      for (k in seq_len(K1)) {
        cols_k <- ((k - 1L) * q + 1L):(k * q)
        add[r, k] <- sum(Znew[r, ] * fit$b[g, cols_k])
      }
    }
  }

  eta <- fhat + add

  exp_eta <- exp(eta)
  exp_eta[is.infinite(exp_eta)] <- 1e50
  denom <- 1 + rowSums(exp_eta)

  prob <- matrix(0, N_new, K)
  prob[, 1:K1] <- exp_eta / denom
  prob[, K] <- 1 / denom
  colnames(prob) <- fit$classes

  pred_idx <- max.col(prob, ties.method = "first")
  fit$classes[pred_idx]
}
