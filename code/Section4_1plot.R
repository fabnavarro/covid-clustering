rm(list=ls())
require(forcats)
library(latex2exp)
library(ggplot2)
load("resultsSection4_1.rda")

###################################################
### Basis Comparison
resultsBasisComparison$bdelta <- as.factor(resultsBasisComparison$bdelta)
resultsBasisComparison$covariateeffect <- as.factor(resultsBasisComparison$covariateeffect)
resultsBasisComparison$model <- as.factor(paste0("scenario ", resultsBasisComparison$model))
resultsBasisComparison$bdelta <- as.factor(resultsBasisComparison$bdelta)
resultsBasisComparison$size <- as.factor(paste0("n-",resultsBasisComparison$size))
resultsBasisComparison$size <- fct_relevel(resultsBasisComparison$size, c("n-50", "n-100", "n-200"))
resultsBasisComparison$basis <- fct_relevel(resultsBasisComparison$basis, c("TI", "polyg", "Bsplines", "Expo"))
levels(resultsBasisComparison$basis) <- c("TI wave", "Poly basis", "Bsplines", "Exp-basis")
resultsBasisComparison$value <- as.numeric(resultsBasisComparison$value)

figure2 <- ggplot(resultsBasisComparison, aes(x=bdelta, y=value))+
  scale_fill_grey(start = 0.4, end = 0.8) +
  geom_boxplot(aes(fill=size))  +
  facet_grid(model ~ basis,scale="free")+
  labs(fill = "") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  theme_grey(base_size = 22) +labs(x = TeX("$\\Delta$"), y="Adjusted Rand Index")+theme(legend.position="top")
print(figure2)

###################################################
### Covariate Comparison
resultsCovariateComparison <- na.omit(resultsCovariateComparison)
resultsCovariateComparison$bdelta <- as.factor(resultsCovariateComparison$bdelta)
resultsCovariateComparison$covariateeffect <- as.factor(resultsCovariateComparison$covariateeffect)
resultsCovariateComparison$model <- as.factor(paste0("scenario ", resultsCovariateComparison$model))
resultsCovariateComparison$bdelta <- as.factor(resultsCovariateComparison$bdelta)
resultsCovariateComparison$size <- as.factor(paste0("n-",resultsCovariateComparison$size))
resultsCovariateComparison$size <- fct_relevel(resultsCovariateComparison$size, c("n-50", "n-100", "n-200"))
resultsCovariateComparison$value <- as.numeric(resultsCovariateComparison$value)
resultsCovariateComparison <- resultsCovariateComparison[which(resultsCovariateComparison$cov!="oracle"),]
resultsCovariateComparison$cov <- as.factor(as.character(resultsCovariateComparison$cov ))
resultsCovariateComparison$cov <- fct_relevel(resultsCovariateComparison$cov, c("cov",  "no-cov"))
levels(resultsCovariateComparison$cov) <- c("Proposed", "No Cov")

figure3 <- ggplot(resultsCovariateComparison, aes(x=cov, y=value))+
  scale_fill_grey(start = 0.4, end = 0.8) +
  geom_boxplot(aes(fill=size))  +
  facet_grid(model ~ covariateeffect,scale="free")+
  labs(fill = "Sample size") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  theme_grey(base_size = 22) +labs(x = "Method", y="Adjusted Rand Index")+
  theme(legend.position="top",axis.text.x = element_text(angle = 65, hjust = 1))
figure3

###################################################
### Clustering Comparison
resultsClusteringComparison <- na.omit(resultsClusteringComparison)
resultsClusteringComparison$bdelta <- as.factor(resultsClusteringComparison$bdelta)
resultsClusteringComparison$covariateeffect <- as.factor(resultsClusteringComparison$covariateeffect)
resultsClusteringComparison$model <- as.factor(paste0("scenario ", resultsClusteringComparison$model))
resultsClusteringComparison$bdelta <- as.factor(resultsClusteringComparison$bdelta)
resultsClusteringComparison$size <- as.factor(paste0("n-",resultsClusteringComparison$size))
resultsClusteringComparison$size <- fct_relevel(resultsClusteringComparison$size, c("n-50", "n-100", "n-200"))
resultsClusteringComparison$value <- as.numeric(resultsClusteringComparison$value)
levels(resultsClusteringComparison$clustering) <- c("HDclassif", "Sparse Nonparametric", "Full Nonparametric")
resultsClusteringComparison$clustering<- fct_relevel(resultsClusteringComparison$clustering, c("Sparse Nonparametric", "HDclassif", "Full Nonparametric"))

figure4 <- ggplot(resultsClusteringComparison, aes(x=size, y=value))+
  scale_fill_grey(start = 0.4, end = 0.8) +
  geom_boxplot(aes(fill=clustering))  +
  facet_grid(.~model ,scale="free")+
  labs(fill = "Mixture model") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  theme_grey(base_size = 22) +labs(x = "Method", y="Adjusted Rand Index")+
  theme(legend.position="top",axis.text.x = element_text(angle = 65, hjust = 1))
figure4

