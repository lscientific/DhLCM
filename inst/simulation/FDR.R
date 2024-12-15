library(foreach)
library(doParallel)

# source the defined functions
source('./R/DhLCM.R')

set.seed(123)
K <- 3
J <- 1000
ratio <- 5
N <- floor(J / ratio)
J_unif <- floor(J/20)

cat("K =", K, "\n")
cat("N =", N, "\n")
cat("J =", J, "\n")
cat("J_unif =", J_unif, "\n")

T0 <- matrix(rbeta(J*K, 0.1, 1), J, K)
for (j in 1:J_unif) {
  T0[j, ] <- rep(runif(1, 0.2, 2/3), K)
}

#### ground truth with degree ####
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
R0[R0 > 1] <- 1

# replication
myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)

set.seed(123)
T1 <- Sys.time()

results <- list()
for(i_level in 1:4) {
  print(i_level)
  level <- c(0.01, 0.05, 0.1, 0.2)[i_level]
  results[[i_level]] <- foreach(rep = 1:500, .packages=c('RcppHungarian', 'rARPACK')) %dopar% {
    R <- matrix(NA, N, J)
    for(i in 1:N) {
      for(j in 1:J) {
        R[i, j] <- rbinom(1, 1, R0[i, j])
      }
    }
    svd_res <- svd(R, nu=K, nv=K)
    U <- svd_res$u
    U_hetero <- heteroPCA(R, K, T0=20)
    
    # with heteroPCA and with L2 norm
    est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
                     S0=S0, clustering_only=F)
    T_hat <- est_res$T_hat
    sigma2_hat <- est_res$sigma2_hat
    
    indices_pos <- which(apply(T_hat, 1, function(x) all(x > 0)))
    len_pos <- length(which(indices_pos %in% 1:J_unif))
    
    T_stat <- c()
    for(j in indices_pos) {
      T_1 <- -1
      for(k1 in 1:(K-1)) {
        for(k2 in (k1+1):K) {
          T_1  <- max(T_1 , (T_hat[j, k1] - T_hat[j, k2])^2 / (sigma2_hat[j, k1] + sigma2_hat[j, k2]))
        }
      }
      T_stat <- c(T_stat, T_1)
    }
    pvalues <- c()
    for (t in T_stat) {
      pvalues <- c(pvalues, 1 - (pchisq(t, df=1))^3)
    }
    fdrs <- p.adjust(pvalues, method="BH")
    fdr <- fdrs[which(indices_pos %in% 1:J_unif)]
    
    false_d <- sum(fdr < level)
    true_d <- sum(fdrs[(J_unif+1):length(fdrs)] < level)
    
    type_I <- sum(fdr < level) / len_pos
    FWER <- sum(pvalues[which(indices_pos %in% 1:J_unif)] < level / len_pos) >= 1
    
    list(false_d, true_d, len_pos, type_I, FWER)
  }
}
T2 <- Sys.time()
print(T2-T1)
stopCluster(myCluster)

#### results
false = true = fdp = len = type_I = FWER <- rep(NA, 4)
for(i_level in 1:4) {
  false_d <- sapply(results[[i_level]], "[[", 1)
  true_d <- sapply(results[[i_level]], "[[", 2)
  false[i_level] <- mean(false_d)
  true[i_level] <- mean(true_d)
  fdp[i_level] <- mean(false_d / (false_d + true_d))
  len[i_level] <- mean(sapply(results[[i_level]], "[[", 3))
  type_I[i_level] <- mean(sapply(results[[i_level]], "[[", 4))
  FWER[i_level] <- mean(sapply(results[[i_level]], "[[", 5))
}

# the values below correspond to Table 2 in the paper
print(false)
print(true)
print(fdp)
print(type_I)

