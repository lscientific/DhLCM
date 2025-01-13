## File Name: svd_hetero_compare.R
## File Version: 0.01
## this file compares SVD and HeteroPCA and corresponds Figure 3 and Figure S.5

library(parallel)
library(ggplot2)
library(ggpubr)
library(latex2exp)

# source the defined functions
source('./R/DhLCM.R')

n_rep <- 100 # number of replications
J <- 3000 # number of items
N <- 500 # number of respondents
K <- 2 # number of latent classes

S0 <- c(rep(1, N/2), rep(2, N/2)) # memberships
theta_1 <- c(rep(0.3, J/2), rep(0.5, J/2))
theta_2 <- c(rep(0.1, J/2), rep(0.06, J/2))
T0 <- cbind(theta_1, theta_2) # ground truth item parameters

#### compare clustering: Figure 3 ####

## when data is Bernoulli-distributed 
f <- function(rep) {
  R0 <- t(matrix(c(rep(theta_1, N/2), rep(theta_2, N/2)), J, N))
  
  Omega0 <- runif(N, 0.5, 1.5)
  W <- rep(NA, K)
  for(k in 1:K) {
    W[k] <- sum(Omega0[S0 == k]^2)
  }
  for(k in 1:K) {
    Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(N/2) # for identifiability
  }
  
  R0 <- diag(Omega0) %*% R0
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rbinom(1, 1, R0[i, j])
    }
  }
  E_mat <- R - R0
  diag_hetero <- diag(E_mat %*% t(E_mat))
  raito_hetero <- max(diag_hetero) / min(diag_hetero)
  
  U0 <- svd(R0, nu=K, nv=K)$u
  U <- svd(R, nu=K, nv=K)$u
  U_hetero <- heteroPCA(R, K, 20)
  
  err_svd <- norm(U %*% t(U) %*% U0 - U0, type = "F")
  err_hetero <- norm(U_hetero %*% t(U_hetero) %*% U0 - U0, type = "F")
  
  res <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Bern', 
               S0=S0, clustering_only=T)
  cls_err_svd <- mean(res$S_hat != S0)
  res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
               S0=S0, clustering_only=T)
  cls_err_hetero <- mean(res$S_hat != S0)
  
  return(c(err_svd, err_hetero, cls_err_svd, cls_err_hetero, raito_hetero))
}

set.seed(1234)
t1 <- Sys.time()
results_bernoulli <- mclapply((1:n_rep), f, mc.cores=7)
print(Sys.time() - t1) # takes about 5 min

## when data is Poisson-distributed 
f <- function(rep) {
  R0 <- t(matrix(c(rep(theta_1, N/2), rep(theta_2, N/2)), J, N))
  
  Omega0 <- runif(N, 0.5, 1.5)
  W <- rep(NA, K)
  for(k in 1:K) {
    W[k] <- sum(Omega0[S0 == k]^2)
  }
  for(k in 1:K) {
    Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(N/2) # for identifiability
  }
  
  R0 <- diag(Omega0) %*% R0
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rpois(1, R0[i, j])
    }
  }
  E_mat <- R - R0
  diag_hetero <- diag(E_mat %*% t(E_mat))
  raito_hetero <- max(diag_hetero) / min(diag_hetero)
  
  U0 <- svd(R0, nu=K, nv=K)$u
  U <- svd(R, nu=K, nv=K)$u
  U_hetero <- heteroPCA(R, K, 20)
  
  err_svd <- norm(U %*% t(U) %*% U0 - U0, type = "F")
  err_hetero <- norm(U_hetero %*% t(U_hetero) %*% U0 - U0, type = "F")
  
  res <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Pois', 
               S0=S0, clustering_only=T)
  cls_err_svd <- mean(res$S_hat != S0)
  res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Pois', 
               S0=S0, clustering_only=T)
  cls_err_hetero <- mean(res$S_hat != S0)
  
  return(c(err_svd, err_hetero, cls_err_svd, cls_err_hetero, raito_hetero))
}

set.seed(123)
t1 <- Sys.time()
results_poisson <- mclapply((1:n_rep), f, mc.cores=7)
print(Sys.time() - t1) # takes about 5 min


#### plot Figure 3 first two columns: clustering 
res_array <- array(NA, dim=c(n_rep, 2))
res_array[, 1] <- sapply(results_poisson, "[[", 1)
res_array[, 2] <- sapply(results_poisson, "[[", 2)
dimnames(res_array) <- list(rep=1:n_rep, method=c('SVD', 'HeteroPCA'))
df_res <- as.data.frame.table(res_array)
p1 <- ggplot(df_res, aes(x = method, y = Freq)) + geom_boxplot() + 
  ylab("Frobenius Norm Error") + xlab(NULL) + ggtitle('Poisson') +
  theme(axis.text=element_text(size=10), plot.title = element_text(size=12))
p1

res_array <- array(NA, dim=c(n_rep, 2))
res_array[, 1] <- sapply(results_poisson, "[[", 3)
res_array[, 2] <- sapply(results_poisson, "[[", 4)
dimnames(res_array) <- list(rep=1:n_rep, method=c('SVD', 'HeteroPCA'))
df_res <- as.data.frame.table(res_array)
p2 <- ggplot(df_res, aes(x = method, y = Freq)) + geom_boxplot() + 
  ylab("Clustering Error") + xlab(NULL) + ggtitle('Poisson') +
  theme(axis.text=element_text(size=10), plot.title = element_text(size=12))
p2

res_array <- array(NA, dim=c(n_rep, 2))
res_array[, 1] <- sapply(results_bernoulli, "[[", 1)
res_array[, 2] <- sapply(results_bernoulli, "[[", 2)
dimnames(res_array) <- list(rep=1:n_rep, method=c('SVD', 'HeteroPCA'))
df_res <- as.data.frame.table(res_array)
p3 <- ggplot(df_res, aes(x = method, y = Freq)) + geom_boxplot() + 
  ylab("Frobenius Norm Error") + xlab(NULL) + ggtitle('Bernoulli') +
  theme(axis.text=element_text(size=10), plot.title = element_text(size=12))
p3 

res_array <- array(NA, dim=c(n_rep, 2))
res_array[, 1] <- sapply(results_bernoulli, "[[", 3)
res_array[, 2] <- sapply(results_bernoulli, "[[", 4)
dimnames(res_array) <- list(rep=1:n_rep, method=c('SVD', 'HeteroPCA'))
df_res <- as.data.frame.table(res_array)
p4 <- ggplot(df_res, aes(x = method, y = Freq)) + geom_boxplot() + 
  ylab("Clustering Error") + xlab(NULL) + ggtitle('Bernoulli') +
  theme(axis.text=element_text(size=10), plot.title = element_text(size=12))
p4


#### compare theta estimation

## when data is Bernoulli-distributed 
f <- function(rep) {
  R0 <- t(matrix(c(rep(theta_1, N/2), rep(theta_2, N/2)), J, N))
  
  Omega0 <- runif(N, 0.5, 1.5)
  W <- rep(NA, K)
  for(k in 1:K) {
    W[k] <- sum(Omega0[S0 == k]^2)
  }
  for(k in 1:K) {
    Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(N/2) # for identifiability
  }
  R0 <- diag(Omega0) %*% R0
  
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rbinom(1, 1, R0[i, j])
    }
  }
  
  U <- svd(R, nu=K, nv=K)$u
  U_hetero <- heteroPCA(R, K, 20)
  
  res_svd <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=F)
  res_hetero <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=F)
  
  err_svd <- max(abs(res_svd$T_hat - T0))
  err_hetero <- max(abs(res_hetero$T_hat - T0))
  
  return(list(err_svd, err_hetero))
}

set.seed(123)
results_theta1 <- mclapply((1:n_rep), f, mc.cores=7)


## when data is Poisson-distributed 
f <- function(rep) {
  R0 <- t(matrix(c(rep(theta_1, N/2), rep(theta_2, N/2)), J, N))
  
  Omega0 <- runif(N, 0.5, 1.5)
  W <- rep(NA, K)
  for(k in 1:K) {
    W[k] <- sum(Omega0[S0 == k]^2)
  }
  for(k in 1:K) {
    Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(N/2) # for identifiability
  }
  R0 <- diag(Omega0) %*% R0
  
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rpois(1, R0[i, j])
    }
  }
  
  U <- svd(R, nu=K, nv=K)$u
  U_hetero <- heteroPCA(R, K, 20)
  
  res_svd <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Pois', 
                   S0=S0, clustering_only=F)
  res_hetero <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Pois', 
                      S0=S0, clustering_only=F)
  
  err_svd <- max(abs(res_svd$T_hat - T0))
  err_hetero <- max(abs(res_hetero$T_hat - T0))
  
  return(list(err_svd, err_hetero))
}

set.seed(123)
results_theta2 <- mclapply((1:n_rep), f, mc.cores=7)


#### plot Figure 3 last column: theta estimation
res_array <- array(NA, dim=c(n_rep, 2))
res_array[, 1] <- sapply(results_theta1, "[[", 1)
res_array[, 2] <- sapply(results_theta1, "[[", 2)
dimnames(res_array) <- list(rep=1:n_rep, method=c('SVD', 'HeteroPCA'))
df_res <- as.data.frame.table(res_array)
p5 <- ggplot(df_res, aes(x = method, y = Freq)) + geom_boxplot() + 
  ylab(TeX("Max Abs Error of $\\Theta$")) + xlab(NULL) + ggtitle('Bernoulli') + 
  theme(axis.text=element_text(size=10), plot.title = element_text(size=12))
p5

res_array <- array(NA, dim=c(n_rep, 2))
res_array[, 1] <- sapply(results_theta2, "[[", 1)
res_array[, 2] <- sapply(results_theta2, "[[", 2)
dimnames(res_array) <- list(rep=1:n_rep, method=c('SVD', 'HeteroPCA'))
df_res <- as.data.frame.table(res_array)
p6 <- ggplot(df_res, aes(x = method, y = Freq)) + geom_boxplot() + 
  ylab(TeX("Max Abs Error of $\\Theta$")) + xlab(NULL) + ggtitle('Poisson') + 
  theme(axis.text=element_text(size=10), plot.title = element_text(size=12))
p6


ggarrange(p3, p4, p5, p1, p2, p6, ncol=3, nrow=2)
# ggsave("inst/extdata/simulation_figures/svd_hetero_compare.png", width=9, height=5)
# this figure corresponds to Figure 3 in the paper


#### compare inference: Figure S.5 ####
set.seed(12345)
T0 <- cbind(theta_1, theta_2)
# we will conduct hypothesis testing on the first two items
T0[1, ] <- rep(0.5, K)
T0[2, ] <- c(0.3, 0.7)

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
max(R0)

myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)

# when data is Bernoulli-distributed 
results_inf_bern <- foreach(rep = 1:500, .packages=c('RcppHungarian', 'rARPACK')) %dopar% {
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
  
  # SVD
  est_res <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=F)
  T_hat <- est_res$T_hat
  sigma2_hat <- est_res$sigma2_hat
  T_hat_scaled <- matrix(NA, J, K)
  for(j in 1:J) {
    for(k in 1:K) {
      T_hat_scaled[j, k] <- (T_hat[j, k] - T0[j, k]) / sqrt(sigma2_hat[j, k])
    }
  }
  
  T_1_svd <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1_svd  <- max(T_1_svd , (T_hat[1, k1] - T_hat[1, k2])^2 / (sigma2_hat[1, k1] + sigma2_hat[1, k2]))
    }
  }
  
  T_2_svd <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_2_svd  <- max(T_2_svd , (T_hat[2, k1] - T_hat[2, k2])^2 / (sigma2_hat[2, k1] + sigma2_hat[2, k2]))
    }
  }
  
  # hetero
  U_hetero <- heteroPCA(R, K, T0=20)
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=F)
  T_hat <- est_res$T_hat
  sigma2_hat <- est_res$sigma2_hat
  T_hat_scaled <- matrix(NA, J, K)
  for(j in 1:J) {
    for(k in 1:K) {
      T_hat_scaled[j, k] <- (T_hat[j, k] - T0[j, k]) / sqrt(sigma2_hat[j, k])
    }
  }
  
  T_1_hetero <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1_hetero  <- max(T_1_hetero , (T_hat[1, k1] - T_hat[1, k2])^2 / (sigma2_hat[1, k1] + sigma2_hat[1, k2]))
    }
  }
  
  T_2_hetero <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_2_hetero  <- max(T_2_hetero , (T_hat[2, k1] - T_hat[2, k2])^2 / (sigma2_hat[2, k1] + sigma2_hat[2, k2]))
    }
  }
  
  list(T_1_svd, T_2_svd, T_1_hetero, T_2_hetero)
}

# when data is Poisson-distributed 
results_inf_pois <- foreach(rep = 1:500, .packages=c('RcppHungarian', 'rARPACK')) %dopar% {
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
  
  # SVD
  est_res <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Pois', 
                   S0=S0, clustering_only=F)
  T_hat <- est_res$T_hat
  sigma2_hat <- est_res$sigma2_hat
  T_hat_scaled <- matrix(NA, J, K)
  for(j in 1:J) {
    for(k in 1:K) {
      T_hat_scaled[j, k] <- (T_hat[j, k] - T0[j, k]) / sqrt(sigma2_hat[j, k])
    }
  }
  
  T_1_svd <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1_svd  <- max(T_1_svd , (T_hat[1, k1] - T_hat[1, k2])^2 / (sigma2_hat[1, k1] + sigma2_hat[1, k2]))
    }
  }
  
  T_2_svd <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_2_svd  <- max(T_2_svd , (T_hat[2, k1] - T_hat[2, k2])^2 / (sigma2_hat[2, k1] + sigma2_hat[2, k2]))
    }
  }
  
  # hetero
  U_hetero <- heteroPCA(R, K, T0=20)
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Pois', 
                   S0=S0, clustering_only=F)
  T_hat <- est_res$T_hat
  sigma2_hat <- est_res$sigma2_hat
  T_hat_scaled <- matrix(NA, J, K)
  for(j in 1:J) {
    for(k in 1:K) {
      T_hat_scaled[j, k] <- (T_hat[j, k] - T0[j, k]) / sqrt(sigma2_hat[j, k])
    }
  }
  
  T_1_hetero <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1_hetero  <- max(T_1_hetero , (T_hat[1, k1] - T_hat[1, k2])^2 / (sigma2_hat[1, k1] + sigma2_hat[1, k2]))
    }
  }
  
  T_2_hetero <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_2_hetero  <- max(T_2_hetero , (T_hat[2, k1] - T_hat[2, k2])^2 / (sigma2_hat[2, k1] + sigma2_hat[2, k2]))
    }
  }
  
  list(T_1_svd, T_2_svd, T_1_hetero, T_2_hetero)
}

stopCluster()


#### plot Figure S.5
# the figures generated below correspond to Figure S.5 in the Supplementary Material

## when data is Bernoulli-distributed 
T_1_svd <- sapply(results_inf_bern, "[[", 1)
T_2_svd <- sapply(results_inf_bern, "[[", 2)
T_1_hetero <- sapply(results_inf_bern, "[[", 3)
T_2_hetero <- sapply(results_inf_bern, "[[", 4)

pdf("inst/extdata/simulation_figures/Bern_svd_typeI.pdf")
pvalues <- 1 - (pchisq(T_1_svd, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Bernoulli w/ SVD', cex.main=1.5, cex.lab=1.5)
dev.off()

pdf("inst/extdata/simulation_figures/Bern_svd_typeII.pdf")
pvalues <- 1 - (pchisq(T_2_svd, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Bernoulli w/ SVD', cex.main=1.5, cex.lab=1.5)
dev.off()

pdf("inst/extdata/simulation_figures/Bern_hetero_typeI.pdf")
pvalues <- 1 - (pchisq(T_1_hetero, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Bernoulli w/ HeteroPCA', cex.main=1.5, cex.lab=1.5)
dev.off()

pdf("inst/extdata/simulation_figures/Bern_hetero_typeII.pdf")
pvalues <- 1 - (pchisq(T_2_hetero, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Bernoulli w/ HeteroPCA', cex.main=1.5, cex.lab=1.5)
dev.off()


## when data is Poisson-distributed 
T_1_svd <- sapply(results_inf_pois, "[[", 1)
T_2_svd <- sapply(results_inf_pois, "[[", 2)
T_1_hetero <- sapply(results_inf_pois, "[[", 3)
T_2_hetero <- sapply(results_inf_pois, "[[", 4)

pdf("inst/extdata/simulation_figures/Pois_svd_typeI.pdf")
pvalues <- 1 - (pchisq(T_1_svd, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Poisson w/ SVD', cex.main=1.5, cex.lab=1.5)
dev.off()

pdf("inst/extdata/simulation_figures/Pois_svd_typeII.pdf")
pvalues <- 1 - (pchisq(T_2_svd, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Poisson w/ SVD', cex.main=1.5, cex.lab=1.5)
dev.off()

pdf("inst/extdata/simulation_figures/Pois_hetero_typeI.pdf")
pvalues <- 1 - (pchisq(T_1_hetero, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Poisson w/ HeteroPCA', cex.main=1.5, cex.lab=1.5)
dev.off()

pdf("inst/extdata/simulation_figures/Pois_hetero_typeII.pdf")
pvalues <- 1 - (pchisq(T_2_hetero, df=1))^(choose(K, 2))
qqunif(pvalues, logscale = F, cex=0.5, xlim=c(0, 1), ylim=c(0, 1), col='black', lcol="black", 
       main='Poisson w/ HeteroPCA', cex.main=1.5, cex.lab=1.5)
dev.off()


