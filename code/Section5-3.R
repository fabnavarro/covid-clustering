rm(list=ls())
source("CoVidAnalysisHeader.R")

Kselect <- 4
map_colors <- RColorBrewer::brewer.pal(Kselect, "Pastel1")

Containment <- do.call(rbind.data.frame,lapply(1:8, function(u) data.frame(index=names(COVIDdata$supplementary)[u], value=rowMeans(COVIDdata$supplementary[[u]], na.rm = T), cluster=resCOVID$clustering[[Kselect]]$zhat)))
Containment$cluster <- as.factor(Containment$cluster)
Containment$index <- as.factor(Containment$index)



Health <- do.call(rbind.data.frame,lapply(15:22, function(u) data.frame(index=names(COVIDdata$supplementary)[u], value=rowMeans(COVIDdata$supplementary[[u]], na.rm = T), cluster=resCOVID$clustering[[Kselect]]$zhat)))
Health$value <- log(Health$value)
Health$cluster <- as.factor(Health$cluster)
Health$index <- as.factor(Health$index)


indexs <- do.call(rbind.data.frame,lapply(23:33, function(u) data.frame(index=names(COVIDdata$supplementary)[u], value=rowMeans(COVIDdata$supplementary[[u]], na.rm = T), cluster=resCOVID$clustering[[Kselect]]$zhat)))
indexs$cluster <- as.factor(indexs$cluster)
indexs$index <- as.factor(indexs$index)
indexs <- indexs[which(indexs$index%in%c("I1","I3","I4")),]
levels(indexs$index)[c(1,5,7)] <- c("Containment health", "Government response", "Stringency")

legsize <- 25
#png("../paper/appli_indexgrl.png",  width = 1000, height = 1000)
pindexs <- ggplot(indexs, aes(x=index, y=value))+
  scale_fill_manual(values = map_colors)+
  geom_boxplot(aes(fill=cluster))  +
  theme(legend.position="top",axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = legsize),
        axis.text = element_text(size = legsize),
        legend.title = element_text(size=legsize),
        legend.text = element_text(size=legsize))
print(pindexs)
#dev.off()

overall <- rbind.data.frame(Containment[which(Containment$index%in%c("C1","C2","C3")),],
                            Health[which(Health$index%in%c("H2","H3")),])
levels(overall$index)[c(1:3,10:11)] <- c("School closing",
                             "Workplace closing",
                             "Cancel public events",
                             "Testing policy",
                             "Contact tracing")
#png("../paper/appli_index2.png",  width = 1000, height = 1000)
poverall <- ggplot(overall, aes(x=index, y=value))+
  scale_fill_manual(values = map_colors)+
  geom_boxplot(aes(fill=cluster))  +
  theme_grey(base_size = 22) + theme(legend.position="top",axis.text.x = element_text(angle = 45, hjust = 1),
                                     axis.title = element_text(size = legsize),
                                     axis.text = element_text(size = legsize),
                                     legend.title = element_text(size=legsize),
                                     legend.text = element_text(size=legsize))
print(poverall)
#dev.off()
