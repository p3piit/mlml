#' Predict class labels from a fitted GMERF model
#'
#' Generates binary predictions for new clustered data using a fitted
#' Generalized Mixed-Effects Random Forest (GMERF) model. Predictions combine
#' the fixed-effects random forest component with cluster-specific random
#' effects estimated during training.
#'
#' For clusters observed in the training data, the corresponding random effects
#' are added to the fixed-effects prediction. For unseen clusters, only the
#' fixed-effects component is used.
#'
#' @param fit A fitted GMERF model object returned by \code{\link{fit_gmerf}}.
#' @param new_df A data.frame containing new observations for prediction.
#' @param thr Numeric. Classification threshold for converting predicted
#'   probabilities into class labels (default: 0.5).
#' @param random_effects Character. Name of the covariate used for the random
#'   slope in the random-effects design.
#' @param id Character. Name of the cluster identifier column in \code{new_df}.
#'
#' @return A numeric vector of predicted class labels (0/1).
#'
#' @details
#' The linear predictor is computed as:
#' \deqn{\eta = \hat f(X) + Z_i \hat b_i}
#' where \eqn{\hat f(X)} is the random forest prediction and \eqn{Z_i \hat b_i}
#' is the random-effects contribution for cluster \eqn{i}, when available.
#'
#' @family gmerf
#'
#' @seealso \code{\link{fit_gmerf}} for model fitting.
#'
#' @references
#' Hajjem, Ahlem, François Bellavance, and Denis Larocque. 2011.
#' “Mixed Effects Regression Trees for Clustered Data.” \emph{Statistics & Probability Letters} 81 (4): 451–59.
#' \doi{10.1016/j.spl.2010.12.003}.
#'
#' Hajjem, Ahlem, Denis Larocque, and François Bellavance. 2017.
#' “Generalized Mixed Effects Regression Trees.” \emph{Statistics & Probability Letters} 126 (July): 114–18.
#' \doi{10.1016/j.spl.2017.02.033}.
#'
#' Pellagatti, Massimo, Chiara Masci, Francesca Ieva, and Anna M. Paganoni. 2021.
#' “Generalized Mixed-Effects Random Forest: A Flexible Approach to Predict University Student Dropout.”
#' \emph{Statistical Analysis and Data Mining: The ASA Data Science Journal} 14 (3): 241–57.
#' \doi{10.1002/sam.11505}.
#'
#' @examples
#' # Simulate clustered data
#' df <- sim_data_gmert(G = 20, n_i = 30, seed = 1)
#'
#' # Fit GMERF
#' fit <- fit_gmerf(df, ntrees = 200, seed = 1, , max_iter_out = 50)
#'
#' # Predict on training data
#' yhat <- predict_gmerf(fit, df)
#' table(yhat, df$y)
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
predict_gmerf <- function(fit,
                          new_df,
                          thr = 0.5,
                          random_effects = "x1",
                          id = "id"
) {
  # fixed part: regression tree predictions using predictors only
  fhat <- as.numeric(predict(fit$forest, data = new_df)$predictions)
  # construct random effects design for new data: [1, x1]
  Znew <- cbind(1, new_df[random_effects])
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
