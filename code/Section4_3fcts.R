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
require(geoR)

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

rCovidDataGRF <- function(n, noiselevel, bdelta, covariateeffect,phi=1){
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
  coord <- matrix(runif(n*2),n,2)
  noises <- replicate(512,grf(n, cov.pars=c(1, phi/n), grid = coord, mean=0, messages = F)$data)
  for (i in 1:n){
    obs[i,] <- obs[i,] + sqrt(noiselevel) * noises[i,]
    obs[i,] <- obs[i,] *  (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
    mu[i] <- (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
  }
  list(obs=obs, z=z, covariates=covariates, bdelta=bdelta, covariateeffect=covariateeffect, mu = mu, model=noiselevel)
}

rsimu4_1GRF <- function(n,
                     bdelta=0,
                     covariateeffect=0,
                     model=1,
                     phi=1){
  delta <- sample(0:bdelta, n, replace=TRUE)
  size <- 512
  z <- sort(sample(1:3, n, replace=TRUE, prob = c(.5,.25,.25)))
  obs <- matrix(NA, n, size) 
  if (model==1){
    effect <- .5
    fbase <- sin(2.5 * pi * (1:(size))/size)
    fbasebis <- sin((2.5-effect) * pi * (1:(size))/size)
    fbase <- fbase * (fbase>0)
    fbasebis <- fbasebis * (fbasebis>0)
    r <- rbind(fbase, (1+effect) * fbase, fbasebis)
    for (i in 1:n) obs[i,] <- c(rep(0, delta[i]), r[z[i], (1:(size-delta[i]))] )
  }else if (model==2){
    x <- seq(0,1,length.out = size)
    f1 <- sin(5*pi*x)+ sin(2*pi*x)
    fg1 <- t(replicate(sum(z==1),f1+rnorm(size)))
    gamma1 <- gamma2 <-  matrix(0, nrow=4, ncol = 4, byrow = T)
    diag(gamma1) <- c(0.672, 0.228, 0.9, 0.34)
    gamma1[1,2] <- -0.134
    gamma1[2,1] <- 0.364
    diag(gamma2) <- diag(gamma1)
    fg2 <- simul.far(m=sum(z==2), n=size, base=base.simul.far(24, 5),
                     d.rho=diag(c(0.45, 0.9, 0.34, 0.45)), alpha=gamma1, cst1=1)
    fg3 <- simul.far(m=sum(z==3), n=size, base=base.simul.far(24, 5), d.rho=diag(c(0.45, 0.9, 0.34, 0.45)),
                     alpha=gamma2, cst1=1)
    fall <- rbind(fg1,fg2$var,fg3$var)
    for (i in 1:n) obs[i,] <- c(rep(0, delta[i]), fall[i, (1:(size-delta[i]))] )
  } 
  covariates <- matrix(rnorm(n * 2), n, 2)
  obsnorm <- obs
  mu <- rep(NA, n)
  coord <- matrix(runif(n*2),n,2)
  noises <- replicate(size,grf(n, cov.pars=c(1, phi/n), grid = coord, mean=0, messages = F)$data)
  
  for (i in 1:n){
    obs[i,] <- obs[i,] + sqrt(0.5) * noises[i,]
    obsnorm[i,] <- obs[i,]
    obs[i,] <- obs[i,] *  (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
    mu[i] <- (covariateeffect * ((sum(covariates[i,]/sqrt(2)))**2) + 1) / (1+covariateeffect)
  }
  list(obs=obs, obsnorm=obsnorm, z=z,  covariates=covariates, bdelta=bdelta, covariateeffect=covariateeffect, mu = mu, model=model)
}