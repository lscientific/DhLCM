## File Name: atac_cellxgene.R
## File Version: 0.01
## this file conducts real data analysis on the single cell dataset
## download the single cell data from figshare https://figshare.com/s/9b4d5964af498d167e85


library(ggplot2)
library(GGally)
library(Seurat)
library(fossil)

# source the defined functions
source('./R/DhLCM.R')

# load the full data
path <- '/Users/lingchen/Downloads/' # path to figshare download
full_data <- readRDS(file.path(paste0(path, 'cell_types.rds')))

# check cell types
# cell type abbreviations:
# endothelial cell: EN  
# smooth muscle cell: SM
# stromal cell: ST
# pericyte: PE          
# leukocyte: LE
cell_types <- full_data$cell_type
print(table(cell_types)) 
# K is the number of extreme profiles
K <- length(unique(cell_types))
print(K)

# full gene expression data
# rows are cells and columns are genes
data <- t(as.matrix(full_data$RNA$counts))
rm(full_data)

# sample at most 500 for each cell type to ensure balanced clusters
# R is the final data matrix to be analyzed
set.seed(123)
R <- matrix(NA, nc=ncol(data))
chosen <- c()
for (k in 1:K) {
  indices <- which(cell_types == unique(cell_types)[k])
  idx <- sample(indices, min(length(indices), 500))
  chosen <- c(chosen, idx)
  R <- rbind(R, data[idx, ])
}
R <- R[-1, ]
genes <- colnames(R)
cell_types <- cell_types[chosen]

# check maximum value
print(max(as.vector(R))) # 212
# check sparsity
print(mean(R != 0)) 

N <- nrow(R)
J <- ncol(R)

# SVD
svd_res <- svds(R, K) 
U <- svd_res$u
V <- svd_res$v
U_hetero <- heteroPCA(R, K, 20)


#### plot ####
U_data <- data.frame(U_hetero[, 1:2], Type=as.factor(cell_types))
p_atac <- ggplot(U_data, aes(x=X1, y=X2, color = Type, alpha=0.8)) + scale_alpha(guide = 'none') + geom_point(size=3) + 
  theme(legend.position = "bottom", legend.key.size = unit(0.6, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=16), axis.text=element_blank(), axis.title=element_blank()) + 
  scale_fill_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883')) +
  scale_colour_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883')) + 
  guides(color = guide_legend(nrow = 2))
p_atac
# this figure corresponds to the second figure in Figure 1 in the paper


U_data <- data.frame(U_hetero[, 1:5], Type=as.factor(cell_types))
ggpairs(U_data, aes(colour = Type, alpha=0.8), legend = 1, columns = 1:5, 
        upper  = list(continuous = "blank"), 
        lower = list(continuous = wrap(ggally_points, size = 0.5))) +
  theme(legend.position = "right", legend.key.size = unit(0.3, 'cm'), 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=9),
        axis.text=element_blank(), axis.line=element_blank()) + scale_alpha(guide = 'none') + 
  scale_fill_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883')) +
  scale_colour_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883'))
# this figure corresponds to Figure S.8 in the supplementary Material


#### Seurat ####
map <- list('endothelial cell'=1, 'smooth muscle cell'=2, 'stromal cell'=3, 'pericyte'=4, 'leukocyte'=5)
cell_types_num <- unlist(map[cell_types])
names(cell_types_num) <- cell_types

# create Seurat object and conduct the pre-processing pipeline
atac <- CreateSeuratObject(counts = t(R))
atac <- NormalizeData(atac, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
all.genes <- rownames(atac)
atac <- ScaleData(atac, features=all.genes)
atac@meta.data$type <- cell_types_num
atac <- SetIdent(atac, value='type')
head(Idents(atac))

# find all markers
atac.markers <- FindAllMarkers(atac)
print(head(atac.markers))

pvalues_seurat <- atac.markers$p_val
# use BH to adjust the p-values
fdrs_seurat <- p.adjust(pvalues_seurat, method="BH") 
genes_seurat_rej_bh <- unique(atac.markers$gene[which(fdrs_seurat <= 0.05)])
print(length(genes_seurat_rej_bh)) # 14277


#### clustering ####
# The results correspond to Table 3 in the paper
S0 <- as.integer(cell_types)
head(S0)

set.seed(1234)

# heteroPCA w/ L2 normalization
res1 <- DhLCM(R, K=5, spectral=U_hetero, norm='L2', dist='Pois', 
              T0=20, nstart=100, S0=S0, clustering_only=F)
cat("Classification error =", mean(res1$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res1$S_hat, S0), "\n") 

# svd w/ L2 normalization
res3 <- DhLCM(R, K=5, spectral=U, norm='L2', dist='Pois', 
              T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res3$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res3$S_hat, S0), "\n") 

# heteroPCA w/ SCORE normalization
res11 <- DhLCM(R, K=5, spectral=U_hetero, norm='SCORE', dist='Pois', 
               T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res11$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res11$S_hat, S0), "\n") 

# svd w/ SCORE normalization
res33 <- DhLCM(R, K=5, spectral=U, norm='SCORE', dist='Pois', 
               T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res33$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res33$S_hat, S0), "\n") 

# heteroPCA + w.o. normalization
res2 <- DhLCM(R, K=5, spectral=U_hetero, norm=NULL, dist='Pois', 
              T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res2$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res2$S_hat, S0), "\n") 

# svd w.o. normalization
res4 <- DhLCM(R, K=5, spectral=U, norm=NULL, dist='Pois', 
              T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res4$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res4$S_hat, S0), "\n") 


#### hypothesis testing ####
# the results correspond to Section 6.3 in the paper
T_hat <- res1$T_hat
sigma2_hat <- res1$sigma2_hat
# indices for positive estimated item parameter locations
indices_pos <- which(apply(T_hat, 1, function(x) all(x > 0)))

indices_rej <- c()
T_stat <- c()
pvalues <- c()
for (idx in indices_pos) {
  T_1 <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1  <- max(T_1 , (T_hat[idx, k1] - T_hat[idx, k2])^2 / (sigma2_hat[idx, k1] + sigma2_hat[idx, k2]))
    }
  }
  pvalues <- c(pvalues, 1 - (pchisq(T_1, df=1))^(choose(K, 2)))
  T_stat <- c(T_stat, T_1)
}

# FDR
fdrs <- p.adjust(pvalues, method="BH")
indices_rej <- which(fdrs < 0.05)

length(indices_pos) # 17300
length(indices_rej) # 12045

genes_rej <- genes[indices_rej]
# number of intersected rejected genes comparing our method and Seurat
length(unique(intersect(genes_rej, genes_seurat_rej_bh))) # 8974






