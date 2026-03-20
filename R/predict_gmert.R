#'
#' Predict binary outcomes from a fitted GMERT model
#'
#' Generates cluster-aware predictions for new data using a fitted
#' generalized mixed effects regression tree (GMERT) model. Predictions
#' combine the fixed-effects regression tree component with estimated
#' cluster-specific random effects when available.
#'
#' @param fit A fitted GMERT model object returned by \code{fit_gmert()}.
#' @param new_df A data.frame containing the predictors and cluster identifiers.
#' @param random_effect Character. Name of the random-effects covariate (default: "x1").
#' @param thr Numeric. Classification threshold for predicted probabilities.
#' @param id Character. Name of the cluster identifier column (default: "id").
#'
#' @return Numeric vector of binary predictions (0 or 1).
#'
#' @seealso [fit_gmert()] for the genral GMERT fitting function.
#'
#' @family gmert
#'
#' @examples
#' \dontrun{
#' data <- sim_data_gmert()
#' split <- split_gmert_data(data)
#' fit <- fit_gmert(split$train, , max_iter_out = 50)
#' preds <- predict_gmert(fit, split$test)
#' table(preds, split$test$y)
#' }
#'
#' @references
#' Hajjem, A., Larocque, D., & Bellavance, F. (2017).
#' \emph{Generalized mixed effects regression trees}.
#' Statistics & Probability Letters, 126, 114--118.
#' \doi{10.1016/j.spl.2017.02.033}
#'
#' @import rpart
#' @import ranger
#' @import MASS
#' @import tidyverse
#' @import lme4
#' @import parallel
#' @import ggplot2
#' @import performance
#'
#' @export
predict_gmert <- function(fit, new_df, random_effect = "x1", thr = 0.5, id = "id"
) {
  # fixed part: regression tree predictions using predictors only
  fhat <- as.numeric(predict(fit$tree, new_df))        # tree fitted on y*; use on new data
  # construct random effects design for new data: [1, x1]
  Znew <- cbind(1, new_df[random_effect])
  # container for random effects contributions
  add <- numeric(nrow(new_df))
  # clusters seen during training (so we have estimated b_i)
  clus_fit <- sort(unique(fit$train_ids))
  # map each new cluster id to its index in clus_fit (NA if unseen)
  map <- match(new_df[[id]], clus_fit)
  seen <- !is.na(map)                                  # TRUE for rows belonging to seen clusters
  if (any(seen)) {
    # add random effect contribution: row-wise sum of [1, x1] * b_i
    add[seen] <- rowSums(Znew[seen, , drop = FALSE] *
                           fit$b[map[seen], ])
  }
  # total prediction = fixed part + random effects contribution
  eta   <- fhat + add
  p     <- 1 / (1 + exp(-eta))

  # return 0/1;
  yhat  <- ifelse(p >= thr, 1, 0)
  yhat
}
