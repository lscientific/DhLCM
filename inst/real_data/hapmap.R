library(ggplot2)                  
library(GGally)
library(fossil)

# source the defined functions
source('./R/DhLCM.R')

# load the full data
path <- '/Users/lingchen/Downloads/' # path to figshare download
data_hapmap <- read.csv(paste0(path, "hapmap.traw"), sep="\t")

R <- t(data_hapmap[, -(1:6)])
data_hapmap <- data_hapmap[, 1:6]

N <- nrow(R)
J <- ncol(R)
cat("N =", N, ", J =", J, "\n")

# check sparsity
print(mean(R != 0)) 

data_eth <- read.csv(paste0(path, 'hapmap_eth.csv'))[, 2]
data_eth_factor <- factor(data_eth, 
                          levels = c("ASW", "CEU", "CHB", "CHD", "GIH", "JPT", 
                                     "LWK", "MEX", "MKK", "TSI", "YRI"), 
                          labels = c("ASW", "CEU", "CHB", "CHD", "GIH", "JPT", 
                                     "LWK", "MEX", "MKK", "TSI", "YRI"))
S0 <- as.numeric(factor(data_eth, labels = 1:11))


# SVD
K <- length(unique(data_eth)) # 11
svd_res <- svds(R, K)
U <- svd_res$u
V <- svd_res$v
U_hetero <- heteroPCA(R, K, 20)


#### plot ####
U_data <- data.frame(U_hetero[, c(4, 6)], Ethnicity=as.factor(data_eth_factor))
p_hapmap <- ggplot(U_data, aes(x=X1, y=X2, color = Ethnicity, alpha=0.8)) + scale_alpha(guide = 'none') + geom_point(size=3) + 
  theme(legend.position = "bottom", legend.key.size = unit(0.6, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=16), 
        axis.text=element_blank(), axis.title=element_blank()) + 
  scale_fill_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd')) +
  scale_colour_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd')) + 
  guides(color = guide_legend(nrow = 2))
p_hapmap
# this figure corresponds to the third figure in Figure 1 in the paper

U_data <- data.frame(U_hetero[, 1:7], Ethnicity=data_eth_factor)
ggpairs(U_data, aes(colour = Ethnicity, alpha=0.8), legend = 1, columns = 1:7, 
        upper  = list(continuous = "blank"), 
        lower = list(continuous = wrap(ggally_points, size = 0.5))) +
  theme(legend.position = "right", legend.key.size = unit(0.3, 'cm'), 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=9),
        axis.text=element_blank(), axis.line=element_blank()) + 
  scale_alpha(guide = 'none') + 
  scale_fill_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd')) +
  scale_colour_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd'))
# this figure corresponds to Figure S.7 in the Supplementary Material


#### clustering ####
# The results correspond to Table 3 in the paper
set.seed(1234)

# heteroPCA w/ L2 normalization
res1 <- DhLCM(R, K=K, spectral=U_hetero, norm='L2', dist='Bern', 
              T0=20, nstart=100, S0=S0, clustering_only=F)
cat("Classification error =", mean(res1$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res1$S_hat, S0), "\n") 

# svd w/ L2 normalization
res3 <- DhLCM(R, K=K, spectral=U, norm='L2', dist='Bern', 
              T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res3$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res3$S_hat, S0), "\n") 

# heteroPCA w/ SCORE normalization
res11 <- DhLCM(R, K=K, spectral=U_hetero, norm='SCORE', dist='Bern', 
              T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res11$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res11$S_hat, S0), "\n") 

# svd w/ SCORE normalization
res33 <- DhLCM(R, K=K, spectral=U, norm='SCORE', dist='Bern', 
              T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res33$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res33$S_hat, S0), "\n") 

# heteroPCA + w.o. normalization
res2 <- DhLCM(R, K=K, spectral=U_hetero, norm=NULL, dist='Bern', 
               T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res2$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res2$S_hat, S0), "\n") 

# svd w.o. normalization
res4 <- DhLCM(R, K=K, spectral=U, norm=NULL, dist='Bern', 
               T0=20, nstart=100, S0=S0, clustering_only=T)
cat("Classification error =", mean(res4$S_hat != S0), "\n") 
cat("Rand index =", rand.index(res4$S_hat, S0), "\n") 


