rm(list=ls())
require(Clustfun)
require(FactoMineR)
load("../COVIDfull.rda")
set.seed(123)


res.pca <- PCA(COVIDdata$covariates[,c(1,2,13)], ncp = Inf)
res.pca2 <- PCA(COVIDdata$covariates[,-c(1:5,13)], ncp = Inf)
cov.pca <- cbind(res.pca$ind$coord[,1:2],res.pca2$ind$coord[,1:3])
resCOVID <- clustfun(COVIDdata$deathrate, 1:10, covariates = cov.pca[,c(1,3,4)])

#save(resCOVID, cov.pca, file="CoVidresultsfull.rda")


