log_prac <- function(x) {
  res <- sapply(x, function(y) {
    log(max(y, 1e-320))
  })
  return(res)
}

loglik <- function(R, T_mat, S_vec) {
  res <- 0
  for(i in 1:N) {
    res <- res + sum(R[i,] * log_prac(T_mat[, S_vec[i]])) + sum((1-R[i,]) * log_prac(1 - T_mat[, S_vec[i]]))
  }
  return(res)
}

update_T <- function(R, Z_mat) {
  denom <- colSums(Z_mat) # K x 1
  T_hat <- t(apply(t(R) %*% Z_mat, 1, function(x) x / denom))
  #T_hat[, which(denom==0)] <- 0
  
  return(T_hat)
}

update_S <- function(R, T_mat) {
  T_mat_0 <- apply(T_mat, 2, log_prac)
  T_mat_1 <- apply(1 - T_mat, 2, log_prac)
  S_hat <- apply(R %*% T_mat_0 + (1-R) %*% T_mat_1, 1, which.max)
  
  return(S_hat)
}

LCM_em <- function(Z_init, R, K, max_iter=100, tol=1e-3) { # Z: Z initialization
  N <- nrow(R)
  J <- ncol(R)
  
  # update theta
  T_hat <- update_T(R, Z_init)
  
  # update Z
  S_hat <- update_S(R, T_hat)
  Z_hat <- matrix(0, N, K)
  for(i in 1:N) {
    Z_hat[i, S_hat[i]] <- 1
  }
  
  # log likelihood
  loglik_curr <- loglik(R, T_hat, S_hat)
  loglik_prev <- -1
  loglik_save <- c(loglik_curr)
  
  iter <- 1
  while((iter <= max_iter) & (abs((loglik_curr - loglik_prev) / loglik_prev) > tol)) {
    iter <- iter + 1
    loglik_prev <- loglik_curr
    
    # update theta
    T_hat <- update_T(R, Z_hat)
    #T_hat[which(T_hat==0)] <- 1e-4
    
    # update Z
    S_prev <- S_hat
    S_hat <- update_S(R, T_hat)
    Z_hat <- matrix(0, N, K)
    for(i in 1:N) {
      Z_hat[i, S_hat[i]] <- 1
    }
    
    # log likelihood
    loglik_curr <- loglik(R, T_hat, S_hat)
    loglik_save  <- c(loglik_save, loglik_curr)
  }
  
  return(list(S_hat=S_hat, Z_hat=Z_hat, iter=iter, loglik_save=loglik_save))
}
