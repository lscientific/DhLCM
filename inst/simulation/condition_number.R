library(stats)
library(combinat)
library(ordinal)
library(RcppHungarian)
library(rARPACK)
library(gap)
library(foreach)
library(doParallel)

# source the defined functions
source('./R/DhLCM.R')


K <- 3
J <- 1000
ratio <- 5
N <- floor(J / ratio)
cat("N =", N, "\n")
cat("J =", J, "\n")

#### Bernoulli ####

set.seed(123)

T0 <- matrix(rbeta(J*K, 0.1, 1), J, K)
T0 <- 2/3 * T0
svd_res <- svd(T0)
cn0 <- svd_res$d[1] / svd_res$d[K]
# change cn in 1, 2.29, 2.975, 3.74, 4.59, 5.62, 6.905
# so that cn_new is 1.5, 2, 2.5, 3, 3.5, 4, 4.5

cn <- 6.905
r <- cn / cn0

T0_new <- T0 %*% diag(c(1, 1, 1/r))
T0_new[T0_new > 1] <- 1
T0_new[1, ] <- rep(0.5, K)
T0_new[2, ] <- 0.1 * c(1, 3, 6)
svd_res1 <- svd(T0_new)
cn_new <- svd_res1$d[1] / svd_res1$d[K]
cat("cn0=", cn0, ", cn_new=", cn_new, '\n')

# R0
S0 <- sample(1:K, N, replace=T) # membership
C0 <- table(S0)
Z0 <- matrix(0, N, K)
for (i in 1:N) {
  Z0[i, S0[i]] <- 1
}

Omega0 <- runif(N, 0.5, 1.5)
W <- rep(NA, K)
for(k in 1:K) {
  W[k] <- sum(Omega0[S0 == k]^2)
}
for(k in 1:K) {
  Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(C0[[k]]) # for identifiability
}

R0 <- diag(Omega0) %*% Z0 %*% t(T0_new)
R0[R0>1] <- 1

idx <- which(T0_new > 0.3)[K+1]
k_chosen <- floor(idx / J) + 1
j_chosen <- idx - J * (k_chosen - 1)

# replication
myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)

results <- foreach(rep = 1:500, .packages=c('RcppHungarian', 'rARPACK')) %dopar% {
  cat("rep =", rep, "\n")
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rbinom(1, 1, R0[i, j])
    }
  }
  svd_res <- svd(R, nu=K, nv=K)
  U <- svd_res$u
  V <- svd_res$v
  Sigma <- svd_res$d[1:K]
  
  # estimation
  est_res <- DhLCM(R, K=K, spectral='heteroPCA', norm='L2', dist='Bern', 
                   S0=S0, clustering_only=F)
  T_hat <- est_res$T_hat
  sigma2_hat <- est_res$sigma2_hat
  class_err <- mean(est_res$S_hat != S0)
  # cat("Classification error =", class_err, "\n")
  
  T_hat_scaled <- matrix(NA, J, K)
  for(j in 1:J) {
    for(k in 1:K) {
      T_hat_scaled[j, k] <- (T_hat[j, k] - T0[j, k]) / sqrt(sigma2_hat[j, k])
    }
  }
  
  # type I error for a single theta
  typeI_single <- abs(T_hat_scaled[j_chosen, k_chosen]) > qnorm(0.975)
  
  # type I error for a row in theta
  T_1 <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1  <- max(T_1 , (T_hat[1, k1] - T_hat[1, k2])^2 / (sigma2_hat[1, k1] + sigma2_hat[1, k2]))
    }
  }
  typeI_multiple <- NA
  
  # power for a row in theta
  T_2 <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_2  <- max(T_2 , (T_hat[2, k1] - T_hat[2, k2])^2 / (sigma2_hat[2, k1] + sigma2_hat[2, k2]))
    }
  }
  power <- NA
  
  list(typeI_single, typeI_multiple, power, class_err, T_1, T_2)
}

T_1 <- sapply(results, "[[", 5); T_2 <- sapply(results, "[[", 6)
T_1 <- T_1[!is.na(T_1)]; T_2 <- T_2[!is.na(T_2)]

# type I error and power
# the values printed below correspond to Table S.1 in the Supplementary Material
beta <- (1 - 0.05)^{1 / choose(K, 2)}
print(mean(T_1 > qchisq(beta, 1)))
print(mean(T_2 > qchisq(beta, 1)))


#### Poisson ####
# Mapping of cn_new and cn:
# 1.5: 1.525
# 2: 2.42
# 2.5: 3.095
# 3: 3.76
# 3.5: 4.416
# 4: 5.08

set.seed(123)
T0 <- matrix(rgamma(J*K, 0.4, 1), J, K)
svd_res <- svd(T0)
cn0 <- svd_res$d[1] / svd_res$d[K]

cn <- 5.08
r <- cn / cn0
T0_new <- T0 %*% diag(c(1, 1, 1/r))
T0_new[1, ] <- rep(0.5, K)
T0_new[2, ] <- 0.1 * c(1, 3, 6)
svd_res1 <- svd(T0_new)
cn_new <- svd_res1$d[1] / svd_res1$d[K]
cat("cn0=", cn0, ", cn_new=", cn_new, '\n')

Omega0 <- runif(N, 0.5, 1.5)
S0 <- sample(1:K, N, replace=T) # membership
C0 <- table(S0)
Z0 <- matrix(0, N, K)
for (i in 1:N) {
  Z0[i, S0[i]] <- 1
}
W <- rep(NA, K)
for(k in 1:K) {
  W[k] <- sum(Omega0[S0 == k]^2)
}
for(k in 1:K) {
  Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(C0[[k]]) # for identifiability
}
R0 <- diag(Omega0) %*% Z0 %*% t(T0_new)

idx <- which(T0_new > 1)[K+1]
k_chosen <- floor(idx / J) + 1
j_chosen <- idx - J * (k_chosen - 1)

# replication
myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)

results <- foreach(rep = 1:500, .packages=c('RcppHungarian', 'rARPACK')) %dopar% {
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rpois(1, R0[i, j])
    }
  }
  svd_res <- svd(R, nu=K, nv=K)
  U <- svd_res$u
  V <- svd_res$v
  
  # estimation
  est_res <- DhLCM(R, K=K, spectral='heteroPCA', norm='L2', dist='Pois', 
                   S0=S0, clustering_only=F)
  T_hat <- est_res$T_hat
  sigma2_hat <- est_res$sigma2_hat
  class_err <- mean(est_res$S_hat != S0)
  cat("Classification error =", class_err, "\n")
  
  T_hat_scaled <- matrix(NA, J, K)
  for(j in 1:J) {
    for(k in 1:K) {
      T_hat_scaled[j, k] <- (T_hat[j, k] - T0[j, k]) / sqrt(sigma2_hat[j, k])
    }
  }
  
  # type I error for a single theta
  typeI_single <- abs(T_hat_scaled[j_chosen, k_chosen]) > qnorm(0.975)
  
  # type I error for a row in theta
  T_1 <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1  <- max(T_1 , (T_hat[1, k1] - T_hat[1, k2])^2 / (sigma2_hat[1, k1] + sigma2_hat[1, k2]))
    }
  }
  typeI_multiple <- NA
  
  # power for a row in theta
  T_2 <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_2  <- max(T_2 , (T_hat[2, k1] - T_hat[2, k2])^2 / (sigma2_hat[2, k1] + sigma2_hat[2, k2]))
    }
  }
  power <- NA
  
  list(typeI_single, typeI_multiple, power, class_err, T_1, T_2)
}

T_1 <- sapply(results, "[[", 5); T_2 <- sapply(results, "[[", 6)
T_1 <- T_1[!is.na(T_1)]; T_2 <- T_2[!is.na(T_2)]

# type I error and power
# the values printed below correspond to Table S.1 in the Supplementary Material
beta <- (1 - 0.05)^{1 / choose(K, 2)}
print(mean(T_1 > qchisq(beta, 1)))
print(mean(T_2 > qchisq(beta, 1)))

