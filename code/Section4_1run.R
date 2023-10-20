rm(list=ls())
source("Section4_1fcts.R")
# Sample size
nvals <- c(50, 100, 200)
# Number of replicats
nbrep <- 100
# Results storage
resultsBasisComparison <- resultsCovariateComparison <- resultsClusteringComparison <- resultsCovariateInvestigate <- data.frame()
#load("resultsSection4_1.rda")
# simulation setup
bdeltavals <-  c(0,20,40,60)
covariateeffectvals <-  c(0,0.1,0.2,0.3)
# nvals <- 100
# nbrep=30

for (bdelta in bdeltavals){
  for (n in nvals){
    set.seed(123)
    all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=bdelta, covariateeffect=0, model=1), simplify = F)
    out <- do.call(rbind.data.frame, mclapply(all.ech, oneBasisComparison, mc.cores = 40))
    resultsBasisComparison <- rbind.data.frame(resultsBasisComparison, out)
    set.seed(123)
    all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=bdelta, covariateeffect=0, model=2), simplify = F)
    out <- do.call(rbind.data.frame, mclapply(all.ech, oneBasisComparison, mc.cores = 40))
    resultsBasisComparison <- rbind.data.frame(resultsBasisComparison, out)
    save(resultsBasisComparison, resultsCovariateComparison, resultsCovariateInvestigate, resultsClusteringComparison, file="resultsSection4_1.rda")
  }
}

for (covariateeffect in covariateeffectvals){
  for (n in nvals){
    set.seed(123)
    all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=0, covariateeffect=covariateeffect, model=1), simplify = F)
    out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovariateComparison, mc.cores = 40))
    resultsCovariateComparison <- rbind.data.frame(resultsCovariateComparison, out)
    set.seed(123)
    all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=0, covariateeffect=covariateeffect, model=2), simplify = F)
    out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovariateComparison, mc.cores = 40))
    resultsCovariateComparison <- rbind.data.frame(resultsCovariateComparison, out)
    save(resultsBasisComparison, resultsCovariateComparison, resultsCovariateInvestigate, resultsClusteringComparison, file="resultsSection4_1.rda")
  }
}

for (covariateeffect in covariateeffectvals){
  for (n in nvals){
    set.seed(123)
    all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=0, covariateeffect=covariateeffect, model=1), simplify = F)
    out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovariateInvestigate, mc.cores = 40))
    resultsCovariateInvestigate <- rbind.data.frame(resultsCovariateInvestigate, out)
    set.seed(123)
    all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=0, covariateeffect=covariateeffect, model=2), simplify = F)
    out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovariateInvestigate, mc.cores = 40))
    resultsCovariateInvestigate <- rbind.data.frame(resultsCovariateInvestigate, out)
    save(resultsBasisComparison, resultsCovariateComparison, resultsCovariateInvestigate, resultsClusteringComparison, file="resultsSection4_1.rda")
  }
}

resultsClusteringComparison <- data.frame()
for (n in nvals){
  set.seed(123)
  all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=bdeltavals[2], covariateeffect=covariateeffectvals[2], model=1), simplify = F)
  out <- do.call(rbind.data.frame, mclapply(all.ech, oneClusteringComparison, mc.cores = 40))
  resultsClusteringComparison <- rbind.data.frame(resultsClusteringComparison,  cbind.data.frame( out))
  set.seed(123)
  all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=bdeltavals[2], covariateeffect=covariateeffectvals[2], model=2), simplify = F)
  out <- do.call(rbind.data.frame, mclapply(all.ech, oneClusteringComparison, mc.cores = 40))
  resultsClusteringComparison <- rbind.data.frame(resultsClusteringComparison, cbind.data.frame(out))
  save(resultsBasisComparison, resultsCovariateComparison, resultsCovariateInvestigate,resultsClusteringComparison, file="resultsSection4_1.rda")
}
