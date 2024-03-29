% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Dh-LCM.R
\name{heteroPCA}
\alias{heteroPCA}
\title{HeteroPCA implementation}
\usage{
heteroPCA(R, K, T0)
}
\arguments{
\item{R}{Numeric matrix. The matrix to perform heteroPCA.}

\item{K}{Positive integer. The number of top eigenvectors to be extracted.}

\item{T0}{Positive integer. The number of iterations.}
}
\value{
Numeric matrix \code{U_hat}
}
\description{
This function implements the HeteroPCA algorithm
}
