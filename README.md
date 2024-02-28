# Degree-heterogeneous Latent Class Analysis

## How to Install
```devtools::install_github("lscientific/DhLCM")```

## Overview
Load the R package
```library(DhLCM)```
- ```heteroPCA```: this function performs the HeteroPCA algorithm proposed in [Zhang et al., 2015](https://arxiv.org/abs/1810.08316).
- ```DhLCM```: this function performs kmeans clustering on the top ```K``` eigenvectors/left singular vectors, and estimates the DhLCM model parameters.
