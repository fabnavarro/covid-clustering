rm(list=ls())
source("CoVidAnalysisHeader.R")


##### Table 1: coefficients of the regression with pvalue (variable selection)
# keep <- list(c(1,3,4),c(1,3),c(1,4),c(3,4))
# Table1 <- sapply(keep, regressionpvalue, contj=resCOVID$wave$contj, covariates=cov.pca)
# rownames(Table1) <- c("Environment.PCA1", "Environment.PCA2", "Health.PCA1", "Health.PCA2", "Health.PCA3",  "pvalue")
# xtable(round(Table1,3),digits = 3)

##### Interpretation of the index
df <- data.frame(index=resCOVID$reg$index, mu=resCOVID$reg$effect)
df <- df[which(quantile(df[,1], 0.025)<=df[,1]),]
df <- df[which(quantile(df[,1], 0.975)>=df[,1]),]
p2 <- ggplot(df, aes(x=index, y=mu)) +
  geom_line() +
  theme_grey(base_size = 22) +
  xlab(unname(TeX(c("$X_i^T \\hat{\\gamma}")))) + ylab(unname(TeX(c("$\\hat{\\mu}(X_i^T \\hat{\\gamma})"))))
p3 <- ggplot(df, aes(x=index)) +
  geom_density() +
  theme_grey(base_size = 22) +
  xlab(unname(TeX(c("$X_i^T \\hat{\\gamma}"))))
#png("../paper/appli_coveffec.png", width = 1000, height = 600)
grid.arrange(p2, p3, ncol=2, nrow = 1)
#dev.off()
p3

res.pcaEnvironment <- PCA(COVIDdata$covariates[,c(1,2,13)], ncp =Inf, graph = FALSE)
res.pcaMedical <- PCA(COVIDdata$covariates[,-c(1:5,13)], ncp = Inf, graph = FALSE)
mat <- cbind(as.numeric(cor(res.pcaEnvironment$ind$coord[,1], COVIDdata$covariates[,c(1,2,13,6:12,14)])),
             as.numeric(sapply(c(1,2,13,6:12,14), function(j) cor.test(COVIDdata$covariates[,j], res.pcaEnvironment$ind$contrib[,1])$p.value)),
             as.numeric(cor(res.pcaMedical$ind$coord[,1], COVIDdata$covariates[,c(1,2,13,6:12,14)])),
             as.numeric(sapply(c(1,2,13,6:12,14), function(j) cor.test(COVIDdata$covariates[,j],res.pcaMedical$ind$coord[,1])$p.value)),
             as.numeric(cor(res.pcaMedical$ind$coord[,2], COVIDdata$covariates[,c(1,2,13,6:12,14)])),
             as.numeric(sapply(c(1,2,13,6:12,14), function(j) cor.test(COVIDdata$covariates[,j],res.pcaMedical$ind$coord[,2])$p.value)))
rownames(mat) <- colnames(COVIDdata$covariates[,c(1,2,13,6:12,14)])
xtable(mat)
