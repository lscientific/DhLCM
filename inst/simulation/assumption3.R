## File Name: assumption3.R
## File Version: 0.01
## this file verifies Assumption 3 and corresponds to Figure S.6


library(stats)
library(combinat)
library(ordinal)
library(RcppHungarian)
library(RSpectra)
library(gap)
library(parallel)
library(ggplot2)
library(ggpubr)
library(latex2exp)

# source the defined functions
source('./R/DhLCM.R')

K <- 3
Js <- c(200, 500, 1000, 2000)
n_rep <- 100
alphas <- c(1, 0.95, 1.05)

#### Bernoulli ####
results <- list()
set.seed(1234)
for (i in 1:4) {
  print(i)
  J <- Js[i]
  N <- J / 5
  
  f <- function(rep) {
    # sample ground truth item parameters
    T0 <- matrix(rbeta(J*K, 0.1, 1), J, K)
    T0 <- 2/3 * T0
    
    # sample degree parameters
    Omega0 <- runif(N, 0.5, 1.5)
    
    # sample membership
    S0 <- sample(1:K, N, replace=T) 
    Z0 <- matrix(0, N, K)
    # one-hot encoding for membership
    for (i in 1:N) Z0[i, S0[i]] <- 1
    
    W <- rep(NA, K)
    for(k in 1:K) W[k] <- sum(Omega0[S0 == k]^2)
    C0 <- table(S0)
    
    # without perturbation
    Omega1 <- Omega0
    for(k in 1:K) {
      Omega1[S0 == k] <- Omega1[S0 == k] / sqrt(W[k]) * sqrt(C0[[k]])
    }
    R0 <- diag(Omega1) %*% Z0 %*% t(T0)
    R0[R0 > 1] <- 1
    R <- matrix(NA, N, J)
    for(i in 1:N) {
      for(j in 1:J) {
        R[i, j] <- rbinom(1, 1, R0[i, j])
      }
    }
    res_hetero <- DhLCM(R, K=K, spectral='heteroPCA', norm='L2', dist='Bern', 
                        S0=S0, clustering_only=F)
    err1 <- max(abs(res_hetero$T_hat - T0))
    
    # with perturbation
    Omega2 <- Omega0
    for(k in 1:K) {
      Omega2[S0 == k] <- Omega2[S0 == k] / sqrt(W[k]) * sqrt(C0[[k]]) * alphas[k]
    }
    R0 <- diag(Omega2) %*% Z0 %*% t(T0)
    R0[R0 > 1] <- 1
    R <- matrix(NA, N, J)
    for(i in 1:N) {
      for(j in 1:J) {
        R[i, j] <- rbinom(1, 1, R0[i, j])
      }
    }
    res_hetero <- DhLCM(R, K=K, spectral='heteroPCA', norm='L2', dist='Bern', 
                        S0=S0, clustering_only=F)
    err2 <- max(abs(res_hetero$T_hat - T0))
    
    return(list(err1, err2))
  }
  
  results[[i]] <- mclapply((1:n_rep), f, mc.cores=7)
}

res_array_bern <- array(NA, dim=c(n_rep, 4, 2))
res_array_bern[, 1, 1] <- sapply(results[[1]], "[[", 1)
res_array_bern[, 2, 1] <- sapply(results[[2]], "[[", 1)
res_array_bern[, 3, 1] <- sapply(results[[3]], "[[", 1)
res_array_bern[, 4, 1] <- sapply(results[[4]], "[[", 1)
res_array_bern[, 1, 2] <- sapply(results[[1]], "[[", 2)
res_array_bern[, 2, 2] <- sapply(results[[2]], "[[", 2)
res_array_bern[, 3, 2] <- sapply(results[[3]], "[[", 2)
res_array_bern[, 4, 2] <- sapply(results[[4]], "[[", 2)
dimnames(res_array_bern) <- list(rep=1:n_rep, J=Js, Perturbation=c('w/o', 'w/'))

df_res <- as.data.frame.table(res_array_bern)
ggplot(df_res, aes(x = J, y = Freq, fill=Perturbation)) + geom_boxplot() + 
  ylab(TeX("Max Abs Error of $\\Theta$")) + xlab(NULL) + xlab('J') +
  scale_fill_manual(values=c("#bae1e8", '#e26b57'))

# ggsave("inst/extdata/simulation_figures/assumption3_bern.png", width = 5, height=3.5)
# this figure corresponds to Figure S.6 in the Supplementary Material
