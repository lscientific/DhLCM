#library(rARPACK)
library(ggplot2)                  
library(GGally)

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
print(mean(R != 0)) # 0.268

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
library(ggplot2)
load("../data/hapmap_svd.RData")

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
ggsave("./figs/hapmap_2d.pdf", width=7, height = 6.5)

U_data <- data.frame(U_hetero[, 1:7], Ethnicity=data_eth_factor)
dim(U_data)
ggpairs(U_data, aes(colour = Ethnicity, alpha=0.8), legend = 1, columns = 1:7, 
        upper  = list(continuous = "blank"), 
        lower = list(continuous = wrap(ggally_points, size = 0.5))) +
  theme(legend.position = "right", legend.key.size = unit(0.3, 'cm'), 
        legend.title = element_text(size=10), 
        legend.text = element_text(size=9),
        axis.text=element_blank(), axis.line=element_blank()) + scale_alpha(guide = 'none') + 
  scale_fill_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd')) +
  scale_colour_manual(values=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd'))
ggsave("./figs/HapMap3_pairs.png", width=6, height=5)


library("car")
library(rgl)
scatter3d(U_hetero[,2], U_hetero[,3], U_hetero[,4], surface=FALSE, groups = data_eth_factor, 
          surface.col=c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd'), 
          axis.scales = FALSE, xlab='x1', ylab='x2', zlab='x3')
legend3d(x=0.4, y=0.36, legend = levels(data_eth_factor), magnify=1,
         col = c('#e26b57','#377EC2', '#F2A778', '#F5CB1D', '#c19883', 'pink', '#a17db4', '#5f5f5f', '#7FBFBB', '#C3C1C0', '#fce3bd'),
         pch = 16)
rgl.snapshot(filename = "./figs/hapmap_3d.png")


#### clustering ####
library(fossil)
source('../functions.R')
load("../data/hapmap_svd.RData")

set.seed(1234)

# heteroPCA w/ L2 normalization
res1 <- main_fcn(U_hetero, R, S0, norm='L2', dist='Bern', nstart=100)
cat("Classification error =", mean(res1$S_hat != S0), "\n") 
adj.rand.index(res1$S_hat, S0)

# svd w/ L2 normalization
res3 <- main_fcn(U, R, S0, norm='L2', dist='Bern', nstart=100)
cat("Classification error =", mean(res3$S_hat != S0), "\n")
adj.rand.index(res3$S_hat, S0)

# heteroPCA w/ SCORE normalization
res11 <- main_fcn(U_hetero, R, S0, norm='SCORE', dist='Bern', nstart=100)
cat("Classification error =", mean(res11$S_hat != S0), "\n") 
adj.rand.index(res11$S_hat, S0)

# svd w/ SCORE normalization
res33 <- main_fcn(U, R, S0, norm='SCORE', dist='Bern', nstart=100)
cat("Classification error =", mean(res33$S_hat != S0), "\n")
adj.rand.index(res33$S_hat, S0)

# heteroPCA + w.o. normalization
res2 <- main_fcn(U_hetero, R, S0, norm=NULL, dist='Bern', nstart=100)
cat("Classification error =", mean(res2$S_hat != S0), "\n") 
adj.rand.index(res2$S_hat, S0)

# svd w.o. normalization
res4 <- main_fcn(U, R, S0, norm=NULL, dist='Bern', nstart=100)
cat("Classification error =", mean(res4$S_hat != S0), "\n")
adj.rand.index(res4$S_hat, S0)


#### hypothesis testing ####
K <- ncol(U)
print(K)

T_hat <- res1$T_hat
sigma2_hat <- res1$sigma2_hat
c_K <- 2 * log(choose(K, 2)) - log(log(choose(K, 2))) - log(pi)

indices_pos <- which(apply(T_hat, 1, function(x) all(x > 0)))
T_stat <- c()
for (idx in indices_pos) {
  T_1 <- -1
  for(k1 in 1:(K-1)) {
    for(k2 in (k1+1):K) {
      T_1  <- max(T_1 , (T_hat[idx, k1] - T_hat[idx, k2])^2 / (sigma2_hat[idx, k1] + sigma2_hat[idx, k2]))
    }
  }
  T_1 <- (T_1 - c_K) / 2
  T_stat <- c(T_stat, T_1)
}
length(T_stat)
length(indices_pos)

pvalues <- pgumbel(T_stat, lower.tail = F)
fdrs <- p.adjust(pvalues, method="BH")
summary(fdrs)
print(J) # 274,128
print(length(indices_pos)) # 262,257
sum(fdrs < 0.05) # 260,591
sum(fdrs == 0)

idx_rej_strong <- sort(fdrs, index.return=T)$ix[1:10]
indices_rej_strong <- indices_pos[idx_rej_strong]
T_hat[indices_rej_strong, ]
data_hapmap[indices_rej_strong,]


#dict <- c("ASW"="1", "CEU"="2", "CHB"="3", "CHD"="4", "GIH"="5", "JPT"="6", "LWK"="7", "MEX"="8", "MKK"="9", "TSI"="10", "YRI"="11")
