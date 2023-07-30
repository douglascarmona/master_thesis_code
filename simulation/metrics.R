library("RSKC") # CER can be estimated using RSKC package

#' Computes misclassification rate
#'
#' Missclasification is a commonly used performance measure in subspace
#' clustering. It allows to compare two partitions with the same number of
#' clusters.
#'
#' As getting exact value of misclassification requires checking all
#' permutations and is therefore intrackable even for modest number of clusters,
#' a heuristic approach is proposed. It is assumed that there are K clusters of
#' maximum M elements. Additional requirement is that classes labels are from
#' range [1, K].
#'
#' @param ind integer vector, containing the cluster labels of each case of a partition 1.
#' @param true.ind 	integer vector, containing the cluster labels of each case of a partition 2.
#' @param M an integer, maximum number of elements in the true class.
#' @param clusters na integer, number of clusters.
#' @references {R. Vidal. Subspace clustering. Signal Processing Magazine, IEEE,
#'   28(2):52-68,2011}
#' @return Misclassification rate.
mcr <- function(ind, true.ind, M = max(table(true.ind)), clusters = length(unique(true.ind))) {
  if (length(ind) != length(true.ind)) {
    stop("Partitions have different lengths")
  }
  forbidden <- NULL
  total <- 0
  nG <- max(ind)
  for (i in M:1) { # different concordance levels
    for (j in 1:nG) { # subspaces numbers (found)
      if (sum(j == forbidden) == 0) { # subspace not yet used
        for (cluster in 1:clusters) { # subspaces numbers (true)
          if (sum(j == ind[true.ind == cluster]) == i) {
            total <- total + i
            forbidden <- c(forbidden, j)
            break
          }
        }
      }
    }
  }
  return(1 - total / length(true.ind))
}

#' Computes sensitivity
#' 
#' Requires outliers identified as 0
#'
#' @param ind integer vector, containing the cluster labels of each case of a partition 1.
#' @param true.ind 	integer vector, containing the cluster labels of each case of a partition 2.
sensitivity <- function(ind, true.ind){
    idx_outliers <- which(ind == 0)
    true_idx_outliers <- which(true.ind == 0)
    return(sum(idx_outliers %in% true_idx_outliers) / length(true_idx_outliers))
}


#' Computes specificity
#' 
#' Requires outliers identified as 0
#'
#' @param ind integer vector, containing the cluster labels of each case of a partition 1.
#' @param true.ind 	integer vector, containing the cluster labels of each case of a partition 2.
specificity <- function(ind, true.ind){
    idx_not_outliers <- which(ind != 0)
    true_idx_not_outliers <- which(true.ind != 0)
    return(sum(idx_not_outliers %in% true_idx_not_outliers) / length(true_idx_not_outliers))
}

#' Computes the mean in the 3rd dimension in an array
#'
#' @param metric array with metric results
array_mean <- function(metric){
    return(apply(metric, c(1,2), mean))
}

#' Computes the median in the 3rd dimension in an array
#'
#' @param metric array with metric results
array_median <- function(metric){
    return(apply(metric, c(1,2), median))
}