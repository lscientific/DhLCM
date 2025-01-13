## File Name: hypothesis_testing.R
## File Version: 0.01
## this file conducts Simulation Study II: hypothesis testing
## Corresponds to Figures 4 and S.3 and Table 1


library(stats)
library(RSpectra)
library(RcppHungarian)
library(foreach)
library(doParallel)
library(gap)
library(ordinal)
library(gap)

# source the defined functions
source('./R/DhLCM.R')


#### when the data are Bernoulli-distributed ####

set.seed(123)
K <- 3 # change K between 3 and 10
J <- 1000 # change J here
cat("K =", K, "\n")

myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)

for(ratio in c(10)) {
  N <- floor(J / ratio)
  cat("N =", N, "\n")
  cat("J =", J, "\n")
  
  T0 <- matrix(rbeta(J*K, 0.1, 1), J, K)
  T0 <- 2/3 * T0
  T0[1, ] <- rep(0.5, K)
  if (K == 3) {
    T0[2, ] <- 0.1 * c(1, 3, 6)
  }
  if (K == 10) {
    T0[2, ] <- 0.06 * (1:10)
  }
  
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
  R0 <- diag(Omega0) %*% Z0 %*% t(T0)
  R0[R0>1] <- 1
  summary(as.vector(R0))
  
  idx <- which(T0 > 0.3)[K+1]
  k_chosen <- floor(idx / J) + 1
  j_chosen <- idx - J * (k_chosen - 1)
  
  # replication
  T1 <- Sys.time()
  results <- foreach(rep = 1:500, .packages=c('RcppHungarian', 'RSpectra')) %dopar% {
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
    U_hetero <- heteroPCA(R, K, T0=20)
    est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
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
  T2 <- Sys.time()
  print(T2-T1)
  save(results, file=paste0("inst/extdata/simulation_data/Bern_J=", J, "_N=", N, "_K=", K, ".RData"))
}

stopCluster(myCluster)


#### Bernoulli data: plot Figure 4 and print Table 1 ####
# the figures generated below correspond to Figure 4 in the paper
# the printed values correspond to Table 1 in the paper

# K <- 10
K <- 3
set.seed(123)
c_K <- 2 * log(choose(K, 2)) - log(log(choose(K, 2))) - log(pi)

qchisq_max <- function(p, K) {
  return(qchisq(p^(choose(K, 2)), df=1))
}

# for (J in c(3000, 5000)) {
for (J in c(500, 1000)) {
  for (ratio in c(10)) {
    N <- floor(J / ratio)
    cat("J =", J, " N =", N, "\n")
    load(paste0("inst/extdata/simulation_data/Bern_J=", J, "_N=", N, "_K=", K, ".RData"))
    T_1 <- sapply(results, "[[", 5); T_2 <- sapply(results, "[[", 6)
    print(mean(is.na(T_2)))
    T_1 <- T_1[!is.na(T_1)]; T_2 <- T_2[!is.na(T_2)]
    
    # type I
    beta <- (1 - 0.05)^{1 / choose(K, 2)}
    if (K == 3) {
      beta <- (1 - 0.05)^{1 / choose(K, 2)}
      print(mean(T_1 > qchisq(beta, 1)))
      print(mean(T_2 > qchisq(beta, 1)))
    }
    if (K == 10) {
      print(mean(T_1 > 2*qgumbel(0.95)+c_K))
      print(mean(T_2 > 2*qgumbel(0.95)+c_K))
    }
    
    
    pdf(paste0("inst/extdata/simulation_figures/typeI_K=", K, "_N=", N, "_J=", J, '.pdf'))
    if (K == 3) {
      pvalues <- 1 - (pchisq(T_1, df=1))^(choose(K, 2))
    }
    if (K == 10) {
      pvalues <- pgumbel((T_1 - c_K) / 2, lower.tail = F)
    }
    qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
           main=paste0("K=", K, ", N=", N, ", J=", J), cex.main=1.5, cex.lab=1.5)
    dev.off()
    
    pdf(paste0("inst/extdata/simulation_figures/typeII_K=", K, "_N=", N, "_J=", J, '.pdf'))
    if (K == 3) {
      pvalues <- 1 - (pchisq(T_2, df=1))^(choose(K, 2))
    }
    if (K == 10) {
      pvalues <- pgumbel((T_2 - c_K) / 2, lower.tail = F)
    }
    qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
           main=paste0("K=", K, ", N=", N, ", J=", J), cex.main=1.5, cex.lab=1.5)
    dev.off()
  }
}


#### when the data are Poisson-distributed ####
set.seed(123)
K <- 3 # change K between 3 and 10
J <- 500 # change J here
cat("K =", K, "\n")

myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)

for(ratio in c(10)) {
  N <- floor(J / ratio)
  cat("N =", N, "\n")
  cat("J =", J, "\n")
  
  T0 <- matrix(rgamma(J*K, 0.5, 1), J, K)
  T0[1, ] <- rep(0.5, K)
  if (K == 3) {
    T0[2, ] <- 0.1 * c(1, 3, 6)
  }
  if (K == 10) {
    T0[2, ] <- 0.06 * (1:10)
  }
  
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
  R0 <- diag(Omega0) %*% Z0 %*% t(T0)
  print(summary(as.vector(R0)))
  
  idx <- which(T0 > 1)[K+1]
  k_chosen <- floor(idx / J) + 1
  j_chosen <- idx - J * (k_chosen - 1)
  print(T0[j_chosen, k_chosen])
  
  # replication
  T1 <- Sys.time()
  results <- foreach(rep = 1:500, .packages=c('RcppHungarian', 'RSpectra')) %dopar% {
    R <- matrix(NA, N, J)
    for(i in 1:N) {
      for(j in 1:J) {
        R[i, j] <- rpois(1, R0[i, j])
      }
    }
    svd_res <- svd(R, nu=K, nv=K)
    U <- svd_res$u
    V <- svd_res$v
    Sigma <- svd_res$d[1:K]
    
    # estimation
    U_hetero <- heteroPCA(R, K, T0=20)
    est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Pois', 
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
  T2 <- Sys.time()
  print(T2-T1)
  save(results, file=paste0("inst/extdata/simulation_data/Pois_J=", J, "_N=", N, "_K=", K, ".RData"))
}

stopCluster(myCluster)


#### Poisson data: plot Figure S.3 ####
# the figures generated below correspond to Figure S.3 
# in the Supplementary Material

# K <- 10 # change K between 3 and 10
K <- 3
set.seed(123)
c_K <- 2 * log(choose(K, 2)) - log(log(choose(K, 2))) - log(pi)

# for (J in c(3000, 5000)) {
for (J in c(500, 1000)) {
  for (ratio in c(10)) {
    N <- floor(J / ratio)
    cat("J =", J, " N =", N, "\n")
    load(paste0("inst/extdata/simulation_data/Pois_J=", J, "_N=", N, "_K=", K, ".RData"))
    T_1 <- sapply(results, "[[", 5); T_2 <- sapply(results, "[[", 6)
    print(mean(is.na(T_2)))
    T_1 <- T_1[!is.na(T_1)]; T_2 <- T_2[!is.na(T_2)]
    
    # type I
    beta <- (1 - 0.05)^{1 / choose(K, 2)}
    if (K == 3) {
      beta <- (1 - 0.05)^{1 / choose(K, 2)}
      print(mean(T_1 > qchisq(beta, 1)))
      print(mean(T_2 > qchisq(beta, 1)))
    }
    if (K == 10) {
      print(mean(T_1 > 2*qgumbel(0.95)+c_K))
      print(mean(T_2 > 2*qgumbel(0.95)+c_K))
    }
    
    pdf(paste0("inst/extdata/simulation_figures/Pois_typeI_K=", K, "_N=", N, "_J=", J, '.pdf'))
    if (K == 3) {
      pvalues <- 1 - (pchisq(T_1, df=1))^(choose(K, 2))
    }
    if (K == 10) {
      pvalues <- pgumbel((T_1 - c_K) / 2, lower.tail = F)
    }
    qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
           main=paste0("K=", K, ", N=", N, ", J=", J), cex.main=1.5, cex.lab=1.5)
    dev.off()
    
    pdf(paste0("inst/extdata/simulation_figures/Pois_typeII_K=", K, "_N=", N, "_J=", J, '.pdf'))
    if (K == 3) {
      pvalues <- 1 - (pchisq(T_2, df=1))^(choose(K, 2))
    }
    if (K == 10) {
      pvalues <- pgumbel((T_2 - c_K) / 2, lower.tail = F)
    }
    qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
           main=paste0("K=", K, ", N=", N, ", J=", J), cex.main=1.5, cex.lab=1.5)
    dev.off()
  }
}
