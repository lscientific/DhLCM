# Degree-heterogeneous Latent Class Analysis

## How to Install
```devtools::install_github("lscientific/DhLCM")```

## Overview
Load the R package by running
```library(DhLCM)```
- ```heteroPCA```: this function performs the HeteroPCA algorithm proposed in [Zhang et al., 2015](https://arxiv.org/abs/1810.08316).
- ```DhLCM```: this function performs kmeans clustering on the top ```K``` eigenvectors/left singular vectors of the data matrix, and estimates the DhLCM model parameters.

## Three real datasets
The three real datasets used in our analysis can be found in [https://figshare.com/s/9b4d5964af498d167e85].
