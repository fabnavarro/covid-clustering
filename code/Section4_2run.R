rm(list=ls())
source("Section4_1fcts.R")
source("Section4_2fcts.R")
load("../COVIDfull.rda")
# Sample size
nvals <- c(50, 100, 200)
# Number of replicats
nbrep <- 100
# Results storage
resultsCovidMethodComparison <- resultsCoviddata <- data.frame()
#load("resultsSection4_1.rda")
# simulation setup
bdeltavals <-  c(0,20,40,60)
covariateeffectvals <-  c(0,0.1,0.2,0.3)
# nvals <- 100
# nbrep=30


for (bdelta in bdeltavals[c(1,4)]){
  for (covariateeffect in covariateeffectvals[c(1,4)]){
    for (n in nvals){
      set.seed(123)
      all.ech <- replicate(nbrep, rCovidData(n, noiselevel = 0.5, bdelta=bdelta, covariateeffect=covariateeffect), simplify = F)
      out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovidMethodComparison, mc.cores = 40))
      resultsCoviddata <- rbind.data.frame(resultsCoviddata, out)
      set.seed(123)
      all.ech <- replicate(nbrep, rCovidData(n, noiselevel = 1, bdelta=bdelta, covariateeffect=covariateeffect), simplify = F)
      out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovidMethodComparison, mc.cores = 40))
      resultsCoviddata <- rbind.data.frame(resultsCoviddata, out)
      save(resultsCovidMethodComparison, resultsCoviddata, file="resultsSection4_2.rda")
    }
  }
}

for (bdelta in bdeltavals[c(1,4)]){
  for (covariateeffect in covariateeffectvals[c(1,4)]){
    for (n in nvals){
      set.seed(123)
      all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=bdelta, covariateeffect=covariateeffect, model=1), simplify = F)
      out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovidMethodComparison, mc.cores = 40))
      resultsCovidMethodComparison <- rbind.data.frame(resultsCovidMethodComparison, out)
      set.seed(123)
      all.ech <- replicate(nbrep, rsimu4_1(n, bdelta=bdelta, covariateeffect=covariateeffect, model=2), simplify = F)
      out <- do.call(rbind.data.frame, mclapply(all.ech, oneCovidMethodComparison, mc.cores = 40))
      resultsCovidMethodComparison <- rbind.data.frame(resultsCovidMethodComparison, out)
      save(resultsCovidMethodComparison, resultsCoviddata, file="resultsSection4_2.rda")
    }
  }
}

