
rm(list = ls())
library(microbenchmark)
setwd("~/study/master_thesis_code/")

library("ktaucenters")
library("ktaucenterscpp")
source("./simulation/data_generation.R")
# Dimension where the observation lives
p <- 10

# Number of clusters
clusters_vec <- c(3, 5, 7, 10)

# Contamination level
contamination_level <- 0

# Number of repetitions
seed <- 10

# Distance between cluster centers
cluster_distance <- 20

sample_size <- c()
dimension <- c()
kmeans_time <- c()
ktaucenters_time <- c()
ktaucenterscpp_time <- c()
kmeans_time_std <- c()
ktaucenters_time_std <- c()
ktaucenterscpp_time_std <- c()

for (clt_idx in 1:length(clusters_vec)) {
      # Simulated data
      clusters <- clusters_vec[clt_idx]
      
      simulated_data <- data_simulation(
          seed = seed,
          p = p,
          clusters = clusters,
          contamination_level = contamination_level,
          cluster_distance = cluster_distance,
          amplitude_factor = 5.0
        )
      X <- simulated_data$X
      real_cluster_loc <- simulated_data$cluster_location
      real_cluster_centers <- simulated_data$centers
      
      # Kmeans
      kmeans_clustering <- microbenchmark(
        kmeans(
          x = X,
          centers = clusters,
          iter.max = 100,
          nstart = 10,
          algorithm = "Lloyd"),
        times = 10L)
      
      # ktaucenters
      ktau_clustering <- microbenchmark(
        ktaucenters::ktaucenters(
          X = X,
          K = clusters,
          tolmin = 1e-04,
          NiterMax = 100,
          nstart = 10,
          startWithKmeans = FALSE,
          startWithROBINPD = FALSE,
          cutoff = 0.999
          ),
        times = 10L)
      
      # ktaucenterscpp
      ktaucpp_clustering <- microbenchmark(
        ktaucenterscpp::ktaucenters(
          x = X,
          centers = clusters,
          nstart = 10,
          use_kmeans =  FALSE,
          use_robin = FALSE,
          max_tol =  1e-04,
          cutoff = 0.999
        ),
        times = 10L)
      
      # Processing time
      sample_size <- c(sample_size, nrow(X))
      kmeans_time <- c(kmeans_time, mean(kmeans_clustering$time))
      ktaucenters_time <-c(ktaucenters_time, mean(ktau_clustering$time))
      ktaucenterscpp_time <-c(ktaucenterscpp_time, mean(ktaucpp_clustering$time))
      kmeans_time_std <- c(kmeans_time_std, sd(kmeans_clustering$time))
      ktaucenters_time_std <-c(ktaucenters_time_std, sd(ktau_clustering$time))
      ktaucenterscpp_time_std <-c(ktaucenterscpp_time_std, sd(ktaucpp_clustering$time))
}

plot(
  sample_size,
  kmeans_time / 1000000000,
  ylim = c(0, 16),
  pch = 16,
  col = 2
)
points(sample_size, ktaucenterscpp_time / 1000000000, pch = 16, col = 3)

processing_time <- as.data.frame(
  cbind(
    sample_size,
    kmeans_time/1000000000,
    ktaucenters_time/1000000000,
    ktaucenterscpp_time/1000000000,
    kmeans_time_std/1000000000,
    ktaucenters_time_std/1000000000,
    ktaucenterscpp_time_std/1000000000
    )
  )
write.csv(processing_time, "processing_time.csv", row.names=FALSE)
