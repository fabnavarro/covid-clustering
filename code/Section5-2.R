rm(list=ls())
source("CoVidAnalysisHeader.R")
require(VarSelLCM)

res <- selectNbcompo(resCOVID$obstocluster, 1:10)
Kselect <- res$model@g

xtable(rbind(1:10, sapply(resCOVID$clustering, function(u) u$logSmoothlike)), digits = 0)

coord <- PCA(resCOVID$obstocluster, graph = F)
df <- data.frame(labels=COVIDdata$infos$ID, dim1=coord$ind$coord[,1], dim2=coord$ind$coord[,2], class=as.factor(resCOVID$clustering[[Kselect]]$zhat))
#png("../paper/appliclusterACPmap1.png", width = 1000, height = 600)
p <- ggplot(data = df, aes(x=dim1, y=dim2, label=labels)) +
  theme(legend.position="none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)) +
  geom_text(aes(color=class),size=6) + xlab("First PCA axis")+ ylab("Second PCA axis")
print(p)
#dev.off()

normcurve <- sweep(COVIDdata$deathrate, 1, resCOVID$reg$effect, "/")
out <- (dist(normcurve, method = "DTW"))
m <- matrix(NA, Kselect, Kselect)
for (k in 1:Kselect){
  for (l in 1:Kselect){
    m[k,l] <- mean(as.matrix(out)[which(resCOVID$clustering[[Kselect]]$zhat==k), which(resCOVID$clustering[[Kselect]]$zhat==l)])
  }
}
xtable(m, digits =0)

########## To create Table 5
xtable(cbind(resCOVID$clustering[[Kselect]]$coeff,
             as.numeric(by(rowSums(normcurve),resCOVID$clustering[[Kselect]]$zhat, mean)),
             as.numeric(by(rowSums(normcurve),resCOVID$clustering[[Kselect]]$zhat, sd)),
             as.numeric(by(rowSums(COVIDdata$deathrate),resCOVID$clustering[[Kselect]]$zhat, mean)),
             as.numeric(by(rowSums(COVIDdata$deathrate),resCOVID$clustering[[Kselect]]$zhat, sd)),
             as.numeric(by(resCOVID$reg$effect,resCOVID$clustering[[Kselect]]$zhat, mean)),
             as.numeric(by(resCOVID$reg$effect,resCOVID$clustering[[Kselect]]$zhat, sd))),2)

table(resCOVID$clustering[[Kselect]]$zhat)
paste0(factor(COVIDdata$infos$NameID[which(resCOVID$clustering[[Kselect]]$zhat==1)],nmax = 100))
factor(COVIDdata$infos$NameID[which(resCOVID$clustering[[Kselect]]$zhat==2)],nmax = 100)
factor(COVIDdata$infos$NameID[which(resCOVID$clustering[[Kselect]]$zhat==3)],nmax = 102)
factor(COVIDdata$infos$NameID[which(resCOVID$clustering[[Kselect]]$zhat==4)],nmax = 100)

#boxplot(rowSums(normcurve)~resCOVID$clustering[[Kselect]]$zhat)


ARI(kmeans(rowMeans(COVIDdata$deathrate)/resCOVID$reg$effect,4)$cluster, resCOVID$clustering[[Kselect]]$zhat)


for (k in 1:Kselect){
  who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
  df <- cbind.data.frame(date=rep(1:512, length(who)),
                         deaths=as.numeric(t(normcurve[who,])),
                         region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
  #png(paste0("../paper/applicluster",k,".png"), width = 200, height = 600)
  p1 <- ggplot(df, aes(x = date, y = deaths, colour=region)) +   geom_line(aes()) +
    ylim(min(normcurve), max(normcurve))+  ylab(label="adjusted daily death rate")+ theme(legend.position="none")
  #dev.off()
}

k <- 1
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(normcurve[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p1 <- ggplot(df, aes(x = date, y = deaths, colour=region)) +  ylim(c(0, 55)) +  geom_line(aes()) +    ylab(label="adjusted daily death rate")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                                                     axis.title = element_text(size = 20),
                                                                                                                                                     axis.text = element_text(size = 20))
k <- 2
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(normcurve[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p2 <- ggplot(df, aes(x = date, y = deaths, colour=region)) + ylim(c(0, 55)) + geom_line(aes()) +  ylab(label=" ")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                         axis.title = element_text(size = 20),
                                                                                                                         axis.text = element_text(size = 20))
k <- 3
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(normcurve[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p3 <- ggplot(df, aes(x = date, y = deaths, colour=region)) + ylim(c(0, 55)) +  geom_line(aes()) +   ylab(label=" ")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                           axis.title = element_text(size = 20),
                                                                                                                           axis.text = element_text(size = 20))
k <- 4
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(normcurve[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p4 <- ggplot(df, aes(x = date, y = deaths, colour=region)) + ylim(c(0, 55)) +   geom_line(aes()) +    ylab(label="")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                            axis.title = element_text(size = 20),
                                                                                                                            axis.text = element_text(size = 20))

#require(gridExtra)
#png(paste0("../paper/appliclusters.png"), width = 1000, height = 600)
grid.arrange(p1,p2,p3,p4, nrow=1, widths=rep(1,4)) +theme(legend.position = "none",
                                                          axis.title.x = element_blank(),
                                                          axis.title.y = element_blank(),
                                                          axis.text.x = element_blank(),
                                                          axis.text.y = element_blank(),
                                                          plot.margin = unit(c(3,-5.5,4,3), "mm"))
#dev.off()


k <- 1
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(COVIDdata$deathrate[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p1 <- ggplot(df, aes(x = date, y = deaths, colour=region)) +  ylim(c(0, 55)) +  geom_line(aes()) +    ylab(label="adjusted daily death rate")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                                                     axis.title = element_text(size = 20),
                                                                                                                                                     axis.text = element_text(size = 20))
k <- 2
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(COVIDdata$deathrate[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p2 <- ggplot(df, aes(x = date, y = deaths, colour=region)) + ylim(c(0, 55)) + geom_line(aes()) +  ylab(label=" ")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                         axis.title = element_text(size = 20),
                                                                                                                         axis.text = element_text(size = 20))
k <- 3
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(COVIDdata$deathrate[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p3 <- ggplot(df, aes(x = date, y = deaths, colour=region)) + ylim(c(0, 55)) +  geom_line(aes()) +   ylab(label=" ")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                           axis.title = element_text(size = 20),
                                                                                                                           axis.text = element_text(size = 20))
k <- 4
who <- which(resCOVID$clustering[[Kselect]]$zhat==k)
df <- cbind.data.frame(date=rep(1:512, length(who)), deaths=as.numeric(t(COVIDdata$deathrate[who,])), region=rep(names(resCOVID$clustering[[Kselect]]$zhat)[who], each=512))
p4 <- ggplot(df, aes(x = date, y = deaths, colour=region)) + ylim(c(0, 55)) +   geom_line(aes()) +    ylab(label="")+ theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1),
                                                                                                                            axis.title = element_text(size = 20),
                                                                                                                            axis.text = element_text(size = 20))

#require(gridExtra)
#png(paste0("../paper/appliclustersUnormalized.png"), width = 1000, height = 600)
grid.arrange(p1,p2,p3,p4, nrow=1, widths=rep(1,4)) +theme(legend.position = "none",
                                                          axis.title.x = element_blank(),
                                                          axis.title.y = element_blank(),
                                                          axis.text.x = element_blank(),
                                                          axis.text.y = element_blank(),
                                                          plot.margin = unit(c(3,-5.5,4,3), "mm"))
#dev.off()
