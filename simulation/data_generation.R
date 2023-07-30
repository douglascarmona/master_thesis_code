library("mvtnorm")
library("MASS")

#' Computes mahalanobis distance for a point to each center
#' @param x desired point.
#' @param centers matrix containing cluster centers in each row.
#' @param sigmas numeric vector with variance for each clusters.
mahalanobis_distances <- function(x, centers, sigmas){
  p <- ncol(centers)
  clusters <- nrow(centers)
  distances <- rep(NA, clusters)
  for (idx in 1:clusters){
    distances[idx] <- mahalanobis(x = x, center = centers[idx, ], cov = ginv(diag(sigmas[idx], p)), inverted = TRUE)
  }
  return(distances)
}

#' Checks whether a given point is a outlier
#' @param x desired point.
#' @param centers matrix containing cluster centers in each row.
#' @param sigmas numeric vector with variance for each clusters.
#' @param alpha_level desired confidence region.
is_outlier <- function(x, centers, sigmas, alpha_level) {
  distances <- mahalanobis_distances(x = x, centers = centers, sigmas = sigmas)
  max_distance <- qchisq(1 - alpha_level, df = ncol(centers))
  return(all(distances > max_distance))
}

#' Data generation for simulations
#' @param seed seed for data generation.
#' @param p dimension where the data lives.
#' @param clusters number of clusters in the data.
#' @param contamination_level outlier contamination level.
#' @param cluster_distance separation between real cluster centers.
#' @param max_theta max theta per cluster.
#' @param amplitude_factor amplitude factor for outlier dispersion.
data_simulation <- function(seed, p, clusters, contamination_level, cluster_distance, max_theta = 100, amplitude_factor = 2.0) {
  set.seed(seed)
  if (clusters %% 2 == 0){
    cluster_centers <- cluster_distance * (-0.5 * clusters + seq(from = 1, to = clusters))
  } else {
    cluster_centers <- cluster_distance * ((-0.5 * (clusters - 1) + seq(from = 1, to = clusters)))
  }
  
  # Vector with theta values for each clusters. Theta could be either max_theta or max_theta/2
  theta <- sample(c(1, 0.5), clusters, replace = TRUE) * max_theta
  
  # Data points for each cluster
  n_per_cluster <- min(c(p, 4))*theta
  n_total <- sum(n_per_cluster)
  cluster_location <- rep(1:clusters, n_per_cluster)

  centers <- matrix(cluster_centers, nrow = clusters, ncol = p)
  sigma <- runif(clusters, min = 1, max = 5)

  data_points <- matrix(NA, nrow = n_total, ncol = p)
  for (cluster in 1:clusters) {
    data_points[which(cluster_location == cluster), ] <- rmvnorm(n = n_per_cluster[cluster], mean = centers[cluster,], sigma = diag(sigma[cluster], p))
  }
  
  # Add outliers
  max_bound <- apply(data_points, 2, max)
  min_bound <- apply(data_points, 2, min)
  n_outliers <- floor(contamination_level * n_total) + 1
  
  p_mean <- 0.5 * (min_bound + max_bound)
  delta <- 0.5 * (max_bound - min_bound)
  a <- p_mean - amplitude_factor * delta
  b <- p_mean + amplitude_factor * delta
   
  # Keep outliers outside the confidence region 0.99
  alpha_level <- 1e-2
  outlier_points <- matrix(NA, nrow = n_outliers, ncol = p)
  for (out_idx in 1:n_outliers) {
    # If the point is in the region 1 - alpha_level from one cluster
    # then it is replaced by another point.
    outlier <-  runif(p, min = a, max = b)
    while(!is_outlier(x = outlier, centers = centers, sigmas = sigma, alpha_level = alpha_level)){
      outlier <- runif(p, min = a, max = b)
    }
    outlier_points[out_idx, ] <- outlier
  }
  
  X <- rbind(data_points, outlier_points)
  # Outliers are identified with cluster index = 0
  cluster_location <- c(cluster_location, rep(0, n_outliers))
  
  return(list(X = X, cluster_location = cluster_location, centers = centers, sigma = sigma))
}
