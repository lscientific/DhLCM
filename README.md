# Degree-heterogeneous Latent Class Analysis

## How to Install
```devtools::install_github("lscientific/DhLCM")```

## Overview
Load the R package by running
```library(DhLCM)```
- ```heteroPCA```: this function performs the HeteroPCA algorithm proposed in [Zhang et al., 2015](https://arxiv.org/abs/1810.08316).
- ```DhLCM```: this function performs kmeans clustering on the top ```K``` eigenvectors/left singular vectors of the data matrix, and estimates the DhLCM model parameters.

## Simulation experiments
- Simulation codes are located in `inst/simulation`
- Each file corresponds to a part of simulation experiments, and can be run individually.
>- `clustering.R` conducts clustering analysis on simulated data and generates Figure 2
>- `hypothesis_testing.R` verifies the proposed method in terms of hypothesis testing and generates Figure 4, Figure S.3, Table 1
>- `svd_hetero_compare.R` compares SVD versus HeteroPCA and generates Figure 3, Figure S.5
>- `MLE.R` compares the proposed method, the marginal maximum likelihood (MML) method and the joint maximum likelihood (JML) method, and generates Figure S.1, Figure S.2
>- `assumption3.R` verifies Assumption 3 and generates Figure S.6
>- `condition_number.R` checks the proposed method's robustness to the condition number and generates Table S.1


## Three real data applications
The three pre-processed real datasets used in our analysis can be downloaded from [https://figshare.com/s/9b4d5964af498d167e85]. The raw data can be accessed from [https://legacy.voteview.com/senate112.htm], [https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3], and [https://cellxgene.cziscience.com/collections/d36ca85c-3e8b-444c-ba3e-a645040c6185]
- The pipeline for real data analysis in the paper is located at `inst/real_data`
>- `senate.R` corresponds to Section 6.1. It generates Figure 1, Figure S.9, Table 3
>- `hapmap.R` corresponds to Section 6.2. It generates Figure 1, Figure S.7, Table 3
>- `atac_celxgene.R` corresponds to Section 6.3. It generates Figure 1, Figure S.8, Table 3.
