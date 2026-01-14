#' Compute the F1-score for the majority or minority class
#'
#' Computes the F1-score from a 2x2 confusion matrix, either for the majority
#' class (default) or for the minority class.
#'
#' @param table A 2x2 confusion matrix (e.g., produced by \code{table(y, yhat)}),
#'   where rows correspond to true classes and columns to predicted classes.
#' @param majority Logical. If \code{TRUE}, compute the F1-score for the majority
#'   class; if \code{FALSE}, compute it for the minority class.
#'
#' @return A numeric scalar giving the F1-score for the selected class.
#'
#' @import MASS
#' @import tidyverse
#' @import performance
#'
#' @export
#'
#' @examples
#' # Example confusion matrix
#' cm <- matrix(c(50, 10,
#'                5,  35), nrow = 2, byrow = TRUE)
#'
#' # F1-score for the majority class
#' f1_fun(cm, majority = TRUE)
#'
#' # F1-score for the minority class
#' f1_fun(cm, majority = FALSE)
f1_fun <- function(table,
                   majority = TRUE) {
  mag0 <- sum(table[1,]) > sum(table[2,])
  if ((majority + mag0) != 1) {
    precision <- table[1,1] / (sum(table[1, ]))
    recall    <- table[1,1] / (sum(table[, 1]))
  } else {
    precision <- table[2,2] / (sum(table[2, ]))
    recall    <- table[2,2] / (sum(table[, 2]))
  }
  F1        <- 2 * precision * recall / (precision + recall)
  F1
}
