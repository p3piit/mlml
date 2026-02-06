#'
#' Train/test split by cluster for GMERT data
#'
#' Splits a clustered dataset into training and test sets while preserving
#' the cluster structure. Each cluster contributes the same proportion of
#' observations to the training and test sets.
#'
#' @param df A data.frame containing clustered observations. Must include
#'   a column named \code{id}.
#' @param train_prop Numeric. Proportion of observations per cluster to use
#'   for training.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list with two elements:
#'   \item{train}{Training data.frame}
#'   \item{test}{Test data.frame}
#'
#' @examples
#' data <- sim_data_gmert()
#' split <- split_gmert_data(data)
#' nrow(split$train)
#' nrow(split$test)
#'
#' @family gmert
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
# ---------------------------------------------------------------
# Train/test split per cluster using the sim_data_gmert() output
# ---------------------------------------------------------------
split_gmert_data <- function(df, 
                             train_prop = 0.7, 
                             seed = 123,
                             split_by_cluster = FALSE) {

  # Set random seed to ensure reproducible splits
  set.seed(seed)

  # Initialize vectors to store row indices for train and test sets
  train_idx <- integer(0)
  test_idx  <- integer(0)

  # Extract unique cluster identifiers
  ids <- unique(df$id)

  if (split_by_cluster){
    # Shuffle cluster IDs
    shuffled_ids <- sample(ids)

    # Determine number of clusters for training set
    n_train_clusters <- floor(length(ids) * train_prop)

    # Select cluster IDs for training and test sets
    train_ids <- shuffled_ids[1:n_train_clusters]
    test_ids  <- shuffled_ids[(n_train_clusters + 1):length(ids)]

    # Get row indices for training and test sets based on cluster IDs
    train_idx <- which(df$id %in% train_ids)
    test_idx  <- which(df$id %in% test_ids)
  } else {
  # Loop over each cluster
  for (g in ids) {

    # Row indices corresponding to the current cluster
    idx_g <- which(df$id == g)

    # Number of training observations for this cluster
    n_train_g <- floor(length(idx_g) * train_prop)

    # Randomly sample training indices within the cluster
    train_g <- sample(idx_g, n_train_g)

    # Remaining indices are assigned to the test set
    test_g  <- setdiff(idx_g, train_g)

    # Append to global train/test index vectors
    train_idx <- c(train_idx, train_g)
    test_idx  <- c(test_idx, test_g)
  }
  }
  # Subset the original data into training and test data.frames
  train_df <- df[train_idx, ]
  test_df  <- df[test_idx, ]

  # Return both datasets as a list
  list(train = train_df, test = test_df)
}
