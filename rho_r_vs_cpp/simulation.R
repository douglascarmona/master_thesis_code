

rm(list = ls())
library(microbenchmark)
setwd("~/study/master_thesis_code/")

rhoOpt <- function(x, cc){
  tmp <- x^2 / 2 / (3.25*cc^2)
  tmp2 <- (1.792 - 0.972 * x^2 / cc^2 + 0.432 * x^4 / cc^4 - 0.052 * x^6 / cc^6 + 0.002 * x^8 / cc^8) / 3.25
  tmp[abs(x) > 2*cc] <- tmp2[abs(x) > 2*cc]
  tmp[abs(x) > 3*cc] <- 1
  tmp
}


sample_size <- c(
  10,
  100,
  500,
  1000,
  5000,
  10000,
  50000,
  100000,
  200000,
  500000,
  800000,
  1000000,
  1200000)
c <- 1

r_time <- c()
cpp_time <- c()
r_time_std <- c()
cpp_time_std <- c()
for (i in 1:length(sample_size)){
  set.seed(10)
  data <- rnorm(sample_size[i])
  
  # R time
  r_rho <- microbenchmark(rhoOpt(data, cc=c), unit = "s")
  
  # C++ time
  cpp_rho <- microbenchmark(ktaucenterscpp::rho_opt(data, c), unit = "s")
  
  r_time <- c(r_time, mean(r_rho$time))
  cpp_time <- c(cpp_time, mean(cpp_rho$time))
  r_time_std <- c(r_time_std, sd(r_rho$time))
  cpp_time_std <- c(cpp_time_std, sd(cpp_rho$time))
  
}

plot(
  sample_size,
  r_time,
  #ylim = c(0, 16),
  pch = 16,
  col = 2
)
points(sample_size, cpp_time, pch = 16, col = 3)

processing_time <- as.data.frame(
  cbind(
    sample_size, 
    r_time, 
    cpp_time,
    r_time_std,
    cpp_time_std
    )
  )
write.csv(processing_time, "processing_time_rho_r_vs_cpp.csv", row.names=FALSE)
