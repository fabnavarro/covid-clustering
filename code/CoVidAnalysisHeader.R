load("CoVidresultsfull.rda")
load("COVIDfull.rda")
Sys.setenv("LANGUAGE"="En")
Sys.setlocale("LC_ALL", "en_GB.UTF-8")
require(ggplot2)
#require(ISPR)
require(regpro)
require(zoo)
require(xtable)
library(gridExtra)
library(latex2exp)
library(FactoMineR)
library(dplyr)
#library(maps)
#library(usmap)
require(MASS)
require(dtw)
require(mixtools)
# Function use to investigate the relevance of the covariates
regressionpvalue <- function(contj, covariates, keep){
  covariates <- scale(covariates)
  tmp <- scale(log(sqrt(t(contj))), scale = F)
  cv <- covariates
  target <- tmp[,1]
  for (j in 2:ncol(tmp)){
    cv <- rbind(cv, covariates)
    target <- c(target, tmp[,j])
  }
  siparam <- rep(0, ncol(covariates))
  drop <- NULL
  if (length(keep)!=ncol(covariates)) drop <- (1:ncol(covariates))[-keep]
  siparam[1:length(keep)] <- as.numeric(single.index(cv[,keep,drop=F], target, h = sd(target) * (length(target)**(-1/5))))
  pvalue <- ispr(list(y=target, w=cv[,c(keep,drop)]), list(gamma2=siparam), band = sd(cv[,keep]%*%siparam[1:length(keep)] )*  (length(target)**(-1/5)))$pvalue
  param <- rep(0, ncol(covariates))
  param[keep] <- siparam[1:length(keep)]
  c(param, pvalue)
}


# Function used to select the number of components
selectNbcompo <- function(x, gvals, M=round(nrow(x)**(1/5)), nbcores=1){
  xred <- x
  for (j in 1:ncol(xred)) {
    xred[,j] <- (xred[,j] - min(xred[,j])) / (max(xred[,j]) - min(xred[,j]))
    bounds <- quantile(xred[,j], (1:M)/M)
    tmp = 0
    for (h in 1:length(bounds)) tmp = tmp + (xred[,j]<=bounds[h])
    xred[,j] <- tmp
  }
  xred <- as.data.frame(xred)
  for (j in 1:ncol(xred)) xred[,j] <- as.factor(xred[,j])
  res <- VarSelCluster(xred, gvals, FALSE, nbcores = nbcores)
  model <- res@model
  proba <- matrix(1, nrow(x), 1)
  if (sum(model@omega)>1){
    proba <- npMSL(scale(x[,which(model@omega==1),drop=F]), model@g, h=nrow(x)**(-1/5), samebw = T, verb=F)$posteriors
  }
  list(probapost = proba,
       partition = apply(proba, 1, which.max),
       partitionbad = res@partitions@zMAP,
       model = model)

}

