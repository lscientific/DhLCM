## File Name: clustering.R
## File Version: 0.01
## this file conducts Simulation Study I: Clustering and corresponds to Figure 2


library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)

# source the defined functions
source('./R/DhLCM.R')

set.seed(123)
K <- 10 # number of latent classes
J <- 1000 # number of items
ratio <- 5 # ratio between J and N
n_rep <- 500 # number of replications

cat("K =", K, "\n")
N <- floor(J / ratio) # number of respondents
cat("N =", N, ", J =", J, "\n")

#### ground truth model has degree heterogeneity ####
myCluster <- makeCluster(detectCores()-1, type = "PSOCK")
registerDoParallel(myCluster)
T1 <- Sys.time()
results1 <- foreach(rep = 1:n_rep, .packages=c('RSpectra', 'RcppHungarian'), .errorhandling="pass") %dopar% {
  
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
  
  # rescale degree parameters for model identifiability
  W <- rep(NA, K)
  for(k in 1:K) {
    W[k] <- sum(Omega0[S0 == k]^2)
  }
  for(k in 1:K) {
    Omega0[S0 == k] <- Omega0[S0 == k] / sqrt(W[k]) * sqrt(C0[[k]]) 
  }
  
  # ground truth response probabilities
  R0 <- diag(Omega0) %*% Z0 %*% t(T0)
  R0[R0 > 1] <- 1
  
  # generate response data matrix
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rbinom(1, 1, R0[i, j])
    }
  }
  
  # calculate SVD and HeteroPCA
  svd_res <- svds(R, k=K, nu=K, nv=K)
  U <- svd_res$u
  U_hetero <- heteroPCA(R, K, T0=20)
  
  # with heteroPCA and without norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm=NULL, dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_no <- mean(est_res$S_hat != S0)
  
  # with heteroPCA and with L1 norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L1', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_l1 <- mean(est_res$S_hat != S0)
  
  # with heteroPCA and with L2 norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_l2 <- mean(est_res$S_hat != S0)
  
  # with heteroPCA and with SCORE norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='SCORE', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_s <- mean(est_res$S_hat != S0)
  
  # with SVD and without norm
  est_res <- DhLCM(R, K=K, spectral=U, norm=NULL, dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_no <- mean(est_res$S_hat != S0)
  
  # with SVD and with L1 norm
  est_res <- DhLCM(R, K=K, spectral=U, norm='L1', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_l1 <- mean(est_res$S_hat != S0)
  
  # with SVD and with L2 norm
  est_res <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_l2 <- mean(est_res$S_hat != S0)
  
  # without SVD and with SCORE norm
  est_res <- DhLCM(R, K=K, spectral=U, norm='SCORE', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_s <- mean(est_res$S_hat != S0)
  
  list(class_err_h_no, class_err_h_l1, class_err_h_l2, class_err_h_s, class_err_s_no, class_err_s_l1, class_err_s_l2, class_err_s_s)
}
T2 <- Sys.time()
print(T2-T1) # about 1 min


#### ground truth does not have degree heterogeneity ####
set.seed(123)

T1 <- Sys.time()
results2 <- foreach(rep = 1:n_rep, .packages=c('RSpectra', 'RcppHungarian'), .errorhandling="pass") %dopar% {
  # sample ground truth item parameters
  T0 <- matrix(rbeta(J*K, 0.1, 1), J, K)
  T0 <- 2/3 * T0
  
  # sample ground truth degree parameters
  Omega0 <- runif(N, 0.1, 1.5)
  
  # sample membership
  S0 <- sample(1:K, N, replace=T) # membership
  C0 <- table(S0)
  # one-hot encoding of membership
  Z0 <- matrix(0, N, K)
  for (i in 1:N) {
    Z0[i, S0[i]] <- 1
  }
  
  # ground truth response probabilities
  R0 <- Z0 %*% t(T0)
  R0[R0 > 1] <- 1
  
  # generate response data matrix
  R <- matrix(NA, N, J)
  for(i in 1:N) {
    for(j in 1:J) {
      R[i, j] <- rbinom(1, 1, R0[i, j])
    }
  }
  
  # calculate SVD and HeteroPCA
  svd_res <- svd(R, nu=K, nv=K)
  U <- svd_res$u
  U_hetero <- heteroPCA(R, K, T0=20)
  
  # with heteroPCA and without norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm=NULL, dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_no <- mean(est_res$S_hat != S0)
  
  # with heteroPCA and with L1 norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L1', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_l1 <- mean(est_res$S_hat != S0)
  
  # with heteroPCA and with L2 norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_l2 <- mean(est_res$S_hat != S0)
  
  # with heteroPCA and with SCORE norm
  est_res <- DhLCM(R, K=K, spectral=U_hetero, norm='SCORE', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_h_s <- mean(est_res$S_hat != S0)
  
  # without SVD and without norm
  est_res <- DhLCM(R, K=K, spectral=U, norm=NULL, dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_no <- mean(est_res$S_hat != S0)
  
  # without SVD and with L1 norm
  est_res <- DhLCM(R, K=K, spectral=U, norm='L1', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_l1 <- mean(est_res$S_hat != S0)
  
  # without SVD and with L2 norm
  est_res <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_l2 <- mean(est_res$S_hat != S0)
  
  # without SVD and with SCORE norm
  est_res <- DhLCM(R, K=K, spectral=U, norm='SCORE', dist='Bern', 
                   S0=S0, clustering_only=T)
  class_err_s_s <- mean(est_res$S_hat != S0)
  
  list(class_err_h_no, class_err_h_l1, class_err_h_l2, class_err_h_s, class_err_s_no, class_err_s_l1, class_err_s_l2, class_err_s_s)
}
T2 <- Sys.time()
print(T2-T1) # takes about 1 min

stopCluster(myCluster)


#### Figure 2 plot ####
err_h_no <- sapply(results1, "[[", 1)
err_h_l1 <- sapply(results1, "[[", 2)
err_h_l2 <- sapply(results1, "[[", 3)
err_h_s <- sapply(results1, "[[", 4)
err_s_no <- sapply(results1, "[[", 5)
err_s_l1 <- sapply(results1, "[[", 6)
err_s_l2 <- sapply(results1, "[[", 7)
err_s_s <- sapply(results1, "[[", 8)

dat <- data.frame(err=c(err_h_no, err_h_l2, err_h_s, err_s_no, err_s_l2, err_s_s), 
                  Normalization=as.factor(rep(c(rep('w/o', n_rep), rep('l2', n_rep), 
                                                rep('SCORE', n_rep)), 2)), 
                  Method=as.factor(c(rep('HeteroPCA', 3*n_rep), rep('SVD', 3*n_rep))))
p1 <- ggplot(dat, aes(y = err, fill = Method, x = Normalization)) + geom_boxplot() + 
  ylab(NULL) + theme(legend.position="bottom") + 
  scale_fill_manual(values=c("#FCE3BD", "#BAE1EB")) +
  theme(legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=9.5),
        legend.text = element_text(size=8), 
        axis.title.x=element_blank()) 
p1 # corresponds to the left panel in Figure 2

err_h_no <- sapply(results2, "[[", 1)
err_h_l1 <- sapply(results2, "[[", 2)
err_h_l2 <- sapply(results2, "[[", 3)
err_h_s <- sapply(results2, "[[", 4)
err_s_no <- sapply(results2, "[[", 5)
err_s_l1 <- sapply(results2, "[[", 6)
err_s_l2 <- sapply(results2, "[[", 7)
err_s_s <- sapply(results2, "[[", 8)

dat2 <- data.frame(err=c(err_h_no, err_h_l2, err_h_s, err_s_no, err_s_l2, err_s_s), 
                   Normalization=as.factor(c(rep(c(rep('w/o', n_rep), rep('l2', n_rep), rep('SCORE', n_rep)), 2))), 
                   Method=as.factor(c(rep('HeteroPCA', 3*n_rep), rep('SVD', 3*n_rep))))
p2 <- ggplot(dat2, aes(y = err, fill = Method, x = Normalization)) + geom_boxplot() + 
  ylab(NULL) + theme(legend.position="bottom") + 
  scale_fill_manual(values=c("#FCE3BD", "#BAE1EB")) +
  theme(legend.key.size = unit(1, 'cm'), 
        legend.title = element_text(size=9.5),
        legend.text = element_text(size=8), 
        axis.title.x=element_blank()) 
p2 # corresponds to the right panel in Figure 2

ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend=T, legend = 'bottom')
# ggsave("inst/extdata/simulation_figures/err_strong.pdf", width=7, height=4)
# this figure corresponds to Figure 2 in the paper


