#' Diagonal deletion
#'
#' This function takes in a matrix, and returns the diagonal-deleted matrix
#'
#' @param X Numeric matrix
#' @return A matrix of with diagonals set to 0
#' @export
H_mat <- function(X) {
  diag(X) <- 0
  return(X)
}




#' HeteroPCA implementation
#'
#' This function implements the HeteroPCA algorithm
#'
#' @param R Numeric matrix. The matrix to perform heteroPCA.
#' @param K Positive integer. The number of top eigenvectors to be extracted.
#' @param T0 Positive integer. The number of iterations.
#' @return Numeric matrix \code{U_hat}
#' @export
heteroPCA <- function(R, K, T0) {
  M <- H_mat(R %*% t(R))
  for (t in 1:T0) {
    svd_res <- RSpectra::svds(M, K)
    M_bar <- svd_res$u %*% diag(svd_res$d[1:K]) %*% t(svd_res$v)
    M <- H_mat(M) + diag(diag(M_bar))
  }
  U_hat <- RSpectra::svds(M, K)$u
  return(U_hat)
}




#' Degree-heterogeneous latent class model estimation
#'
#' This function performs k-means clustering on the top \code{K} eigenvectors/left singular vectors, and estimates the DhLCM model parameters
#'
#' @param R Numeric matrix. Data matrix.
#' @param K Positive integer. The number of top eigenvectors/left singular vectors to be extracted.
#' @param spectral Character. One of \code{"heteroPCA"} and \code{"SVD"}. 
#'        Specifies the method to be used to obtain the top \code{K} eigenvectors/left singular vectors. 
#'        \code{"heteroPCA"} implements the heteroPCA method. 
#'        \code{"SVD"} performs ordinary singular vector decomposition.
#' @param norm Character. One of \code{"L2"}, \code{"L1"}, \code{"SCORE"}, and \code{"NA"}. 
#'        Specifies the method to be used for normalization on the eigenvectors/left singular vectors. 
#'        \code{"L2"} performs L2 normalization. 
#'        \code{"L1"} performs L1 normalization.
#'        \code{"SCORE"} performs SCORE normalization.
#'        \code{"NA"} does not perform normalization.
#' @param dist Character. One of \code{"Bernoulli"}, \code{"Binomial"}, and \code{"Poisson"}. 
#'        Specifies the distribution of the ground truth model.
#'        \code{"Bernoulli"} assumes the Bernoulli distribution.
#'        \code{"Binomial"} assumes the Binomial distribution.
#'        \code{"Poisson"} assumes the Poisson distribution.
#' @param T0 Positive integer. The number of iterations for heteroPCA.
#' @param nstart Positive integer. The number of initial starts in the kmeans function.
#' @return Named list. The list is made of:
#' \itemize{
#' \item \code{T_hat} --- Numeric matrix. Estimation of the \eqn{\Theta}{Theta} matrix.
#' \item \code{S_hat} --- Numeric vector. Clustered membership for each subject.
#' \item \code{Z_hat} --- Numeric matrix. Clustered membership for each subject in binary matrix form.
#' \item \code{sigma2_hat} --- Numeric vector (>=0). Asymptotic variance for each element of \code{T_hat}.
#' }
#' @export
DhLCM <- function(R, K, spectral='heteroPCA', norm='L2', dist='Bernoulli', T0=20, nstart=10) { 
  N <- nrow(R)
  J <- ncol(R)
  U <- heteroPCA(R, K, T0)
  
  # k-means clustering
  if (norm == 'NA') {
    kmeans_res <- stats::kmeans(U, centers=K, iter.max=100, nstart=nstart)
  } 
  if (norm == 'L2') {
    U_bar <- t(apply(U, 1, function(u) u / sqrt(sum(u^2))))
    kmeans_res <- stats::kmeans(U_bar, centers=K, iter.max=100, nstart=nstart)
  } 
  if (norm == 'L1') {
    U_bar <- t(apply(U, 1, function(u) u / sum(abs(u))))
    kmeans_res <- stats::kmeans(U_bar, centers=K, iter.max=100, nstart=nstart)
  } 
  if (norm == 'SCORE') {
    U_bar <- t(apply(U, 1, function(u) u / u[1]))
    U_bar <- U_bar[, 2:ncol(U_bar)]
    if (N <= 300) {
      t <- 2 * log(N)
    } else {
      t <- log(N)
    }
    U_bar <- ifelse(abs(U_bar) < t, U_bar, t)
    kmeans_res <- stats::kmeans(U_bar, centers=K, iter.max=100, nstart=nstart)
  }
  S_hat <- kmeans_res$cluster
  centers <- kmeans_res$centers
  
  # parameter estimation
  T_hat = sigma2_hat = Z_hat <- NULL
  Z_hat <- matrix(0, N, K)
  for (i in 1:N) {
    Z_hat[i, S_hat[i]] <- 1
  }
  C_hat <- table(S_hat)
  Omega_hat <- apply(U, 1, function(u) sqrt(sum(u^2))) * sapply(S_hat, function(s) sqrt(C_hat[[s]]))
  T_hat <- t(solve(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% diag(1 / Omega_hat) %*% R)
  T_hat <- ifelse(T_hat > 1, 1, T_hat)
  T_hat <- ifelse(T_hat < 0, 0, T_hat)
  
  # asymptotic variance 
  if ((dist == "Bernoulli") | (dist == "Binomial")) {
    sigma2_hat <- matrix(NA, J, K)
    for(j in 1:J) {
      for(k in 1:K) {
        w_k <- Omega_hat[S_hat == k]
        sigma2_hat[j, k] <- T_hat[j, k] * sum((1 - w_k * T_hat[j, k]) / w_k) / C_hat[[k]]^2 
      }
    }
  }
  if (dist == "Poisson") {
    sigma2_hat <- matrix(NA, J, K)
    for(j in 1:J) {
      for(k in 1:K) {
        w_k <- Omega_hat[S_hat == k]
        sigma2_hat[j, k] <- T_hat[j, k] * sum(1 / w_k)/ C_hat[[k]]^2
      }
    }
  }
  
  return(list(T_hat=T_hat, S_hat=S_hat, Z_hat=Z_hat, sigma2_hat=sigma2_hat))
}

