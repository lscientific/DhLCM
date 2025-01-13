## File Name: MLE.R
## File Version: 0.01
## this file compares JML, MML and the proposed method. 
## Corresponds to Figures S.1 and S.2

library(foreach)
library(doParallel)
library(stats)
library(RcppHungarian)
library(RSpectra)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(poLCA)

# source the defined functions
source('./R/DhLCM.R')
source('./R/em.R')

set.seed(123)
K <- 2
J <- 500
N <- 100
cat("K =", K, "\n")
cat("N =", N, "\n")
cat("J =", J, "\n")

# sample ground truth item parameters
T0 <- matrix(rbeta(J*K, 0.1, 1), J, K)
T0 <- 2/3 * T0

# sample ground truth degree parameters
Omega0 <- runif(N, 0.1, 1.5)

# sample membership
S0 <- sample(1:K, N, replace=T) 
C0 <- table(S0)
# one-hot encoding of membership
Z0 <- matrix(0, N, K)
for (i in 1:N) {
  Z0[i, S0[i]] <- 1
}

# number of replications
n_rep <- 500


#### ground truth model does not have degree heterogeneity ####
# ground truth response probabilities
R0 <- Z0 %*% t(T0)
R0[R0 > 1] <- 1

myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)
T1 <- Sys.time()
results <- foreach(rep = 1:n_rep, .packages=c('RSpectra', 'poLCA')) %dopar% {
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rbinom(1, 1, R0[i, j])
    }
  }
  svd_res <- svd(R, nu=K, nv=K)
  U <- svd_res$u
  
  #### the proposed method ####
  t1 <- Sys.time()
  est_res <- DhLCM(R, K=K, spectral='heteroPCA', norm='L2', dist='Bern', 
                   S0=S0, clustering_only=F)
  class_err_spectral <- mean(est_res$S_hat != S0)
  t2 <- Sys.time()
  t_spectral <- t2 - t1
  print(class_err_spectral)
  
  #### JML and MML ####
  t1 <- Sys.time()
  Z_init <- matrix(0, N, K)
  for (i in 1:N) {
    Z_init[i, sample(1:K, 1)] <- 1
  }
  res_em <- LCM_em(Z_init, R, K, max_iter=100, tol=1e-10)
  t2 <- Sys.time()
  t_jml <- t2 - t1
  
  S_hat <- res_em$S_hat
  class_err_jml <- mean(S_hat != S0)
  class_err_jml <- min(class_err_jml, 1 - class_err_jml)

  #### mml
  t1 <- Sys.time()
  data <- data.frame(1.0 + R)
  f <- as.formula(paste("cbind(", paste0(colnames(data), collapse = ","),")~ 1"))
  res_mml <- poLCA(f, data=data, nclass=K)
  t2 <- Sys.time()
  t_mml <- t2 - t1
  
  S_hat <- res_mml$predclass
  class_err_mml <- mean(S_hat != S0)
  class_err_mml <- min(class_err_mml, 1-class_err_mml)
  
  list(class_err_spectral, class_err_jml, class_err_mml, t_spectral, t_jml, t_mml)
}
T2 <- Sys.time()
print(T2-T1)

stopCluster(myCluster)



#### Figure S.1 plot ####
err_spectral <- sapply(results, "[[", 1)
err_jml <- sapply(results, "[[", 2)
err_mml <- sapply(results, "[[", 3)
t_spectral <- sapply(results, "[[", 4)
t_jml <- sapply(results, "[[", 5)
t_mml <- sapply(results, "[[", 6)
n_rep <- length(err_spectral)

dat <- data.frame(err=c(err_spectral, err_jml, err_mml), 
                  Method=factor(c(rep('Proposed', n_rep), rep('JML', n_rep), rep('MML', n_rep)), levels=c('MML', 'JML', 'Proposed')))
p1 <- ggplot(dat, aes(y = err, x = Method)) + geom_boxplot() + ggtitle('Clustering Error') + 
  ylab(NULL) + theme(legend.position="bottom") + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=9.5),
        legend.text = element_text(size=8), 
        axis.title.x=element_blank()) 
p1

dat <- data.frame(time=c(t_spectral, t_jml, t_mml), 
                  Method=factor(c(rep('Proposed', n_rep), rep('JML', n_rep), rep('MML', n_rep)), levels=c('MML', 'JML', 'Proposed')))
p2 <- ggplot(dat, aes(y = time, x = Method)) + geom_boxplot() + 
  ylab(NULL) + theme(legend.position="bottom") + ggtitle('Computation Time') + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=9.5),
        legend.text = element_text(size=8), 
        axis.title.x=element_blank()) 
p2

ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend=T, legend = 'bottom')
ggsave("inst/extdata/simulation_data/err_mle_no_degree.pdf", width=7, height=3.5)
# this figure corresponds to Figure S.1 in the Supplementary Material


#### ground truth with degree ####
W <- rep(NA, K)
for(k in 1:K) {
  W[k] <- sum(Omega0[S0 == k]^2)
}
for(k in 1:K) {
  Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(C0[[k]]) # for identifiability
}
R0 <- diag(Omega0) %*% Z0 %*% t(T0)
R0[R0 > 1] <- 1


myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)
T1 <- Sys.time()
results <- foreach(rep = 1:n_rep, .packages=c('rARPACK', 'poLCA')) %dopar% {
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rbinom(1, 1, R0[i, j])
    }
  }
  svd_res <- svd(R, nu=K, nv=K)
  U <- svd_res$u
  
  #### sprectral ####
  t1 <- Sys.time()
  est_res <- DhLCM(R, K=K, spectral='heteroPCA', norm='L2', dist='Bern', 
                   S0=S0, clustering_only=F)
  class_err_spectral <- mean(est_res$S_hat != S0)
  t2 <- Sys.time()
  t_spectral <- t2 - t1
  print(class_err_spectral)
  
  #### JML and MML ####
  t1 <- Sys.time()
  Z_init <- matrix(0, N, K)
  for (i in 1:N) {
    Z_init[i, sample(1:K, 1)] <- 1
  }
  res_em <- LCM_em(Z_init, R, K, max_iter=100, tol=1e-10)
  t2 <- Sys.time()
  t_jml <- t2 - t1
  
  S_hat <- res_em$S_hat
  class_err_jml <- mean(S_hat != S0)
  class_err_jml <- min(class_err_jml, 1 - class_err_jml)
  
  #### mml
  t1 <- Sys.time()
  data <- data.frame(1.0 + R)
  f <- as.formula(paste("cbind(", paste0(colnames(data), collapse = ","),")~ 1"))
  res_mml <- poLCA(f, data=data, nclass=K)
  t2 <- Sys.time()
  t_mml <- t2 - t1
  
  S_hat <- res_mml$predclass
  class_err_mml <- mean(S_hat != S0)
  class_err_mml <- min(class_err_mml, 1-class_err_mml)

  list(class_err_spectral, class_err_jml, class_err_mml, t_spectral, t_jml, t_mml)
}
T2 <- Sys.time()
print(T2-T1)

stopCluster(myCluster)


#### Figure S.2 plot ####
err_spectral <- sapply(results, "[[", 1)
err_jml <- sapply(results, "[[", 2)
err_mml <- sapply(results, "[[", 3)
t_spectral <- sapply(results, "[[", 4)
t_jml <- sapply(results, "[[", 5)
t_mml <- sapply(results, "[[", 6)
n_rep <- length(err_spectral)

dat <- data.frame(err=c(err_spectral, err_jml, err_mml), 
                  Method=factor(c(rep('Proposed', n_rep), rep('JML', n_rep), rep('MML', n_rep)), levels=c('MML', 'JML', 'Proposed')))
p1 <- ggplot(dat, aes(y = err, x = Method)) + geom_boxplot() + ggtitle('Clustering Error') + 
  ylab(NULL) + theme(legend.position="bottom") + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=9.5),
        legend.text = element_text(size=8), 
        axis.title.x=element_blank()) 
p1

dat <- data.frame(time=c(t_spectral, t_jml, t_mml), 
                  Method=factor(c(rep('Proposed', n_rep), rep('JML', n_rep), rep('MML', n_rep)), levels=c('MML', 'JML', 'Proposed')))
p2 <- ggplot(dat, aes(y = time, x = Method)) + geom_boxplot() + 
  ylab(NULL) + theme(legend.position="bottom") + ggtitle('Computation Time') + 
  theme(legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=9.5),
        legend.text = element_text(size=8), 
        axis.title.x=element_blank()) 
p2

ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend=T, legend = 'bottom')
ggsave("inst/extdata/simulation_data/err_mle_w_degree.pdf", width=7, height=3.5)
# this figure corresponds to Figure S.2 in the Supplementary Material


