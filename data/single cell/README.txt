
Single cell sequencing data
Reference: https://www.sciencedirect.com/science/article/pii/S2211124722017302
Data access: https://cellxgene.cziscience.com/collections/d36ca85c-3e8b-444c-ba3e-a645040c6185
Sample size: 18,315 cells with 19,298 genes

========================

Demonstration in R:

library(Seurat)

data <- readRDS(file.path('./cell_types.rds'))
cell_types <- data$cell_type # cell types of length 18,315
counts <- t(as.matrix(data$RNA$counts)) # count data with dim 18,315 x 19,298