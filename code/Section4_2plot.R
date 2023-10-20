rm(list=ls())
require(forcats)
library(latex2exp)
library(ggplot2)
load("resultsSection4_2.rda")

###################################################
resultsCovidMethodComparison <- na.omit(resultsCovidMethodComparison)
resultsCovidMethodComparison$bdelta <- as.factor(paste0("Delta-",resultsCovidMethodComparison$bdelta))
resultsCovidMethodComparison$covariateeffect <- as.factor(paste0("q-",resultsCovidMethodComparison$covariateeffect))
resultsCovidMethodComparison$model <- as.factor(paste0("scenario ", resultsCovidMethodComparison$model))
resultsCovidMethodComparison$size <- as.factor(paste0("n-",resultsCovidMethodComparison$size))
resultsCovidMethodComparison$size <- fct_relevel(resultsCovidMethodComparison$size, c("n-50", "n-100", "n-200"))
levels(resultsCovidMethodComparison$basis)[6] <- "Proposed method"
resultsCovidMethodComparison$basis <- fct_relevel(resultsCovidMethodComparison$basis, c("Proposed method", "NMF", "DTW-spectral", "DTW-ward", "DTW-complete", "DTW-single"))
resultsCovidMethodComparison$value <- as.numeric(resultsCovidMethodComparison$value)

#png("simu4_2Comparison.png", width = 1000, height = 600)
figure5 <- ggplot(resultsCovidMethodComparison, aes(x=size, y=value))+
  scale_fill_grey(start = 0.4, end = 0.8) +
  geom_boxplot(aes(fill=basis))  +
  facet_grid(model ~ bdelta + covariateeffect,scale="free")+
  labs(fill = "Method") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  theme_grey(base_size = 22) +labs(x ="Sample size", y="Adjusted Rand Index")+
  theme(legend.position="top")
print(figure5)
#dev.off()
#
# ###################################################
resultsCoviddata$bdelta <- as.factor(paste0("Delta-",resultsCoviddata$bdelta))
resultsCoviddata$covariateeffect <- as.factor(paste0("q-",resultsCoviddata$covariateeffect))
resultsCoviddata$model <- as.factor(paste0("noise variance ", resultsCoviddata$model))
resultsCoviddata$size <- as.factor(paste0("n-",resultsCoviddata$size))
resultsCoviddata$size <- fct_relevel(resultsCoviddata$size, c("n-50", "n-100", "n-200"))
resultsCoviddata$value <- as.numeric(resultsCoviddata$value)
levels(resultsCoviddata$basis)[6] <- "Proposed method"
resultsCoviddata$basis <- fct_relevel(resultsCoviddata$basis, c("Proposed method", "NMF", "DTW-spectral", "DTW-ward", "DTW-complete", "DTW-single"))

#png("simu4_2Covid.png", width = 1000, height = 600)
figure6 <- ggplot(resultsCoviddata, aes(x=size, y=value))+
  scale_fill_grey(start = 0.4, end = 0.8) +
  geom_boxplot(aes(fill=basis))  +
  facet_grid(model ~ bdelta + covariateeffect,scale="free")+
  labs(fill = "Method") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  theme_grey(base_size = 22) +labs(x ="Sample size", y="Adjusted Rand Index")+
  theme(legend.position="top")
figure6
#dev.off()

