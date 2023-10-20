# Data are built from the files COVID-19_LUT.csv and Policy.rds downloaded from https://github.com/CSSEGISandData/COVID-19_Unified-Dataset
# The file UEcountry contains the ISO code of the UE countries + UK
rm(list=ls())
require(zoo)
 giveData <- function(selectedregions, dateStart, dateEnd, Policy, Covariates){
  # keep the informations for the considered regions
  Policy <- Policy[which(Policy$ID%in%selectedregions),]
  Policy$ID <- as.character(Policy$ID)
  Policy <- Policy[order(Policy$ID),]
  Covariates <- Covariates[which(Covariates$ID%in%selectedregions),]
#  Covariates <- Covariates[which(duplicated(Covariates[,1])==F),]
  Covariates$ID <- as.character(Covariates$ID)
  Covariates <- Covariates[order(Covariates$ID),]
  # check that we have the information for all the considered regions for all the dates
  dates <- as.character(unique(Policy$Date))
  nbdates <- length(dates)
  COVIDdata <- list()
  COVIDdata$supplementary <- list()
  vbles <- unique(Policy$PolicyType)
  for (j in 1:length(vbles)){
    tmp <- Policy[which(Policy$PolicyType==vbles[j]),]
    tmp <- tmp[order(tmp$ID),]
    COVIDdata$supplementary[[j]] <- matrix(as.numeric(tmp$PolicyValue), length(unique(tmp$ID)), nbdates, byrow = TRUE,
                             dimnames = list(as.character(unique(tmp$ID)),as.character(unique(tmp$Date))))
    # keep the records from the period of interest
    COVIDdata$supplementary[[j]] <- COVIDdata$supplementary [[j]][,which(dates==dateStart):which(dates==dateEnd)]
  }
  names(COVIDdata$supplementary) <- vbles
  print(vbles)
  COVIDdata$deathrate <- COVIDdata$supplementary$Deaths
  for (j in 2:ncol(COVIDdata$deathrate)) COVIDdata$deathrate[,j] <-  COVIDdata$supplementary$Deaths[,j] -  COVIDdata$supplementary$Deaths[,j-1]

  COVIDdata$deathrate <- as.matrix(sweep(COVIDdata$deathrate, 1, infos$Population, "/")) * (10**6)
  # # Moving average on 7 days
  COVIDdata$deathrate <- t(apply(COVIDdata$deathrate, 1, rollmean, k=7))
  COVIDdata$covariates <- Covariates[,c("PM2.5_PopWtd", "NO2_PopWtd", "WorldPop_Density", "Diabetes", "Obesity", "Smoking", "COPD", "CVD",  "HIV", "Hypertension", "WorldPop_65")]
  COVIDdata$covariates <- scale(COVIDdata$covariates)
  COVIDdata
}


##################
# UEISO1_3C for UE countries + UK
UEISO1_3C <- c("AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GRC", "HUN", "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROU", "SVK",
               "SVN", "ESP", "SWE", "GBR")
# Code for US states
USAID <- c(paste0("US0", 1:9),paste0("US", 10:90))
infos <- read.csv("COVID-19_LUT.csv")
infos <- infos[c(which((infos$ISO1_3C %in%UEISO1_3C + (infos$Level=="Country")) == 2), which(infos$ID%in%USAID)), ]
# drop: American Samoa, Guam, Northern Mariana Islands, Puerto Rico, Virgin Islands
infos <- infos[-c(80:84),]
# print the set of the considered regions
require(stringr)
infos$NameID <- as.character(infos$NameID)
infos$NameID <- sapply(str_split(infos$NameID, ", "), function(u) u[[1]])
infos$ID <- as.character(infos$ID)
infos <- infos[order(infos$ID),]
infos$Regions <- infos$NameID
infos$Regions[which(infos$Level!="State")] <- as.character(infos$ISO1_3C)[which(infos$Level!="State")]
# Curves
Policy <- readRDS("Policy.rds")
# Covariates
Covariates <- readRDS("COVID-19_Static.rds")

##################
COVIDdata <- giveData(infos$ID, "2020-02-27", "2021-07-28", Policy, Covariates)
COVIDdata$infos <- infos
colnames(COVIDdata$deathrate)
#save(COVIDdata, file="../builddata-Aout2021/COVID.rda")
