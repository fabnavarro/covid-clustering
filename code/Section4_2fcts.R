require(Clustfun)
require(fda)
require(mixtools)
require(HDclassif)
require(parallel)
require(icamix)
require(VGAM)
library(rwavelet)
library(far)
library(dtw)
library(NMF)
library(cluster)
library(sClust)

oneCovidMethodComparison <- function(ech){
  zTInp <- clustfun(ech$obs, 3, type = "TI", covariates = ech$covariates)$clustering[[1]]$zhat
  distdtw <- dist(ech$obs,method="DTW")
  outward <- agnes(distdtw, metric = "ward")
  zdtwward <- cutree(outward,3)
  outsingle <- agnes(distdtw, metric = "single")
  zdtwsingle <- cutree(outsingle,3)
  outcomplete <- agnes(distdtw, metric = "complete")
  zdtwcomplete <- cutree(outcomplete,3)
  sim <- 1-as.matrix(distdtw)/max(as.matrix(distdtw))
  outdtwspec <- spectralPAM(sim, 3)
  if (any(ech$obs<0)) ech$obs <- ech$obs - min(ech$obs) + 0.00001
  best <- 0
  for (k in 1:10){
    a <- nmf(ech$obs, k)
    zMNF <- kmeans(a@fit@W, 3)$cluster
    tmp <- ARI(ech$z, zMNF)  
    if (tmp>best) best <- tmp
  }
  # output
  data.frame(bdelta=rep(ech$bdelta,6),
             model=rep(ech$model, 6),
             covariateeffect=rep(ech$covariateeffect, 6),
             size=rep(nrow(ech$obs), 6),
             basis=c("Prop", "DTW-ward", "DTW-single", "DTW-complete", "DTW-spectral", "NMF"),
             value=c(ARI(ech$z, zTInp), ARI(ech$z, zdtwward), ARI(ech$z, zdtwsingle), ARI(ech$z, zdtwcomplete), ARI(ech$z, outdtwspec$cluster), best)
  )
}

rCovidData <- function(n, noiselevel, bdelta, covariateeffect){
  load("COVIDfull.rda")
  centers <- COVIDdata$deathrate[c(13,14,18),]
  prop <- rep(1/3,3)
  z <- sample(1:3, n, replace=TRUE, prob=prop)
  covariates <- matrix(rnorm(n * 2), n, 2)
  obs <- centers[z,]
  for (i in 1:n){
    delta <- sample(0:bdelta,1)
    obs[i,] <- c(rep(0,delta), obs[i,1:(ncol(obs)-delta)])
  }
  mu <- rep(NA, n)
  for (i in 1:n){
    obs[i,] <- obs[i,] + sqrt(noiselevel) * rnorm(ncol(obs))
    obs[i,] <- obs[i,] *  (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
    mu[i] <- (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
  }
  list(obs=obs, z=z, covariates=covariates, bdelta=bdelta, covariateeffect=covariateeffect, mu = mu, model=noiselevel)
}
