## File Name: senate.R
## File Version: 0.01
## this file conducts real data analysis on the senate dataset
## download the senate data from figshare https://figshare.com/s/9b4d5964af498d167e85


library(ggplot2)                  
library(GGally)
library(fossil)

# source the defined functions
source('./R/DhLCM.R')

# load the full data
path <- '/Users/lingchen/Downloads/' # path to figshare download
data <- read.csv(paste0(path, "sen112kh.csv"))

# remove the first row: the president's votes
data <- data[-1, -1]
party <- data$PARTY.CODE
R <- data[, grep("roll", colnames(data) ) ]

# check sparsity
print(mean(R != 0)) 

# check missingness
print(apply(R, 1, function(x) any(x == -9)))  # almost every senator has missing votes
print(sum(apply(R, 2, function(x) any(x == -9))))  # every vote has missing values

# remove parties that are not Dem. or Repub.
idx_rmv <- which(!party %in% c('100', '200'))
# remove senators with missing rate larger than 0.1
idx_rmv <- c(idx_rmv, which(apply(R, 1, function(x) mean(x==-9)) > 0.1))
R <- R[-idx_rmv, ] 

party <- as.character(party[-idx_rmv])
dict <- c('200'='Repub.', '100'='Dem.')
party <- dict[party]
print(length(party))
print(dim(R))

# imputation
set.seed(123)
prob <- apply(R, 1, function(x) { sum(x==1) / sum(x!=-9) })
for(i in 1:nrow(R)) {
  R[i, which(R[i, ] == -9)] <- rbinom(sum(R[i, ] == -9), 1, prob[i])
}
R <- as.matrix(R)

N <- nrow(R)
J <- ncol(R)
cat("N =", N, ", J =", J, "\n")

# SVD
K <- 2
svd_res = svds(R, K)
U <- svd_res$u
V <- svd_res$v
U_hetero <- heteroPCA(R, K, 20)

#### plot ####
U_data <- data.frame(U_hetero)
U_data$Party <- factor(party)
p_senate <- ggplot(U_data, aes(x=X1, y=X2, color = Party, alpha=0.8)) + scale_alpha(guide = 'none') + geom_point(size=3) + 
  theme(legend.position = "bottom", legend.key.size = unit(0.6, 'cm'), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=16), axis.text=element_blank(), 
        axis.title=element_blank()) + 
  scale_colour_manual(values=c('#e26b57','#377EC2'))
p_senate
# this figure corresponds to the first figure in Figure 1 in the paper


#### clustering ####
# The results correspond to Table 3 in the paper
S0 <- ifelse(party == 'Dem.', 1, 2)

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


#### hypothesis testing ####
# the results correspond to Section 6.1 in the paper
T_hat <- res1$T_hat
sigma2_hat <- res1$sigma2_hat
Z_hat <- res1$Z_hat

# indices for positive estimated item parameter locations
indices_pos <- which(apply(T_hat, 1, function(x) all(x > 0)))

T_stat <- c()
for (idx in indices_pos) {
  T_stat <- c(T_stat, (T_hat[idx, 1] - T_hat[idx, 2])^2 / (sigma2_hat[idx, 1] + sigma2_hat[idx, 2]))
}

pvalues <- pchisq(T_stat, df=1, lower.tail = F)
fdrs <- p.adjust(pvalues, method="BH")

indices_rej <- indices_pos[which(fdrs <= 0.05)]
indices_no_rej <- setdiff(indices_pos, indices_rej)

print(length(indices_pos))
print(length(indices_rej)) 
print(length(indices_no_rej))

desc <- read.csv(paste0(path, "s112desc.csv"))
idx_rej_strong <- sort(fdrs, index.return=T)$ix[1:10]
T_stat[idx_rej_strong]
T_hat[indices_pos[idx_rej_strong],]
fdrs[idx_rej_strong]



#### heatmap ####
desc <- desc[indices_pos,]

T_hat <- res1$T_hat[indices_pos, ]
K <- ncol(T_hat)
head(T_hat)
T_hat_scaled <- t(apply(T_hat, 1, function(x) x / sum(abs(x))))
T_diff <- apply(T_hat_scaled, 1, function(x) sum(abs(x - rep(1/K, K))))
items <- sort(T_stat, index.return=T, decreasing=T)$ix[1:20]
pvalues <- pchisq(T_stat, df=1, lower.tail = F)
fdrs <- p.adjust(pvalues, method="BH")

hm <- heatmap(T_hat[items, ])
items_cluster <- hm$rowInd
sort(hm$colInd, index.return=T)$ix

#names <- rownames(T_hat)[items[items_cluster]]
#names <- sapply(names, function(x) substr(x, 5, nchar(x)))
items_coded <- factor(rownames(T_hat)[items[items_cluster]], levels=rownames(T_hat)[items[items_cluster]])

T_data <- data.frame(value=round(as.vector(T_hat[items[items_cluster],]), 3), Item=items_coded, 
                     Profile=as.factor(rep(c('Democrat', 'Republican'), each=length(items))))

ggplot(T_data, aes(x = Item, y = Profile, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = '#F2A778') +
  scale_y_discrete(expand=c(0, 0)) + scale_x_discrete(expand=c(0, 0)) + ylab("") + xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10), 
        panel.border = element_blank(), panel.grid = element_blank(),
        panel.spacing = element_blank(), line = element_blank(), 
        panel.background = element_blank()) 
# this figure corresponds to Figure S.9 in the supplementary Material

