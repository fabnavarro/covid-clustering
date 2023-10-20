rm(list=ls())
load("CoVidresultsfull.rda")
load("COVIDfull.rda")
Sys.setenv("LANGUAGE"="En")
Sys.setlocale("LC_ALL", "en_GB.UTF-8")
require(ggplot2)
require(zoo)
##################################################################################
# Curves description

data <- data.frame(
  day = as.Date("2020-03-14") + rep(0:511,6),
  Region=as.factor(rep(c("Austria","Italy", "New Hampshire", "Pennsylvania","Costa Rica","Brazil"), each=512)),
  value=c(as.numeric(COVIDdata$deathrate[1,]),
          as.numeric(COVIDdata$deathrate[27,]),
          as.numeric(COVIDdata$deathrate[72,]),
          as.numeric(COVIDdata$deathrate[81,]),
          as.numeric(COVIDdata$deathrate[9,]),
          as.numeric(COVIDdata$deathrate[4,])))
##### To create Figure 1
p <- ggplot(data, aes(x=day, y=value)) +
  theme_grey(base_size = 22) +
  theme(legend.position="top")+
  geom_line(aes(color = Region, linetype = Region)) +
  xlab("") + ylab("Daily CoVid-19 death rate per million")+scale_x_date(date_labels = "%Y %b %d")
print(p)
