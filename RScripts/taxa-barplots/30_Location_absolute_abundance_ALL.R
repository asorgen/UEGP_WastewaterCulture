#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-08-21
#Description: 

rm(list=ls())


## .) Set libraries
library(ggplot2)
# library(grid)
library(scales)
library(grDevices)



## .) Set directories
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="PlotPrep", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")




## .) Create abundance plot
taxa_input=c("Filtered_Family")


dFrame=read.table(paste0(inputPath, taxa_input, '_Abundance_by_Location_Antibiotic_Figs.tsv'),sep='\t',header = TRUE)

dFrame2 <- dFrame[dFrame$Percentage >= 0.01,]
taxa <- unique(dFrame2$Class_Family)
dFrame <- dFrame[dFrame$Class_Family %in% taxa,]

colorIDs=read.table(paste0(inputPath,'colorIDs.tsv'),sep='\t',header = TRUE)

taxaPresent <- unique(dFrame$Class_Family)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]
palette <- colors$hexID

dFrame$Location <- factor(dFrame$Location, levels = c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))

plot <- ggplot(dFrame, aes(fill=Class_Family, y=Average, x=Location)) +
  geom_bar(stat="identity") +
  labs(x="Location", y="Average Abundance (CFU/mL)")+
  scale_fill_manual(values=palette)+
  theme_bw()+
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 18))+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 22))+
  theme(axis.text.y = element_text(size = 22))+
  scale_x_discrete(labels =c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))+
  guides(fill=guide_legend(title="Class - Family",ncol = 1))+
  scale_y_continuous(labels = scientific)

ggsave(plot, file=paste0(output, taxa_input, "_absAbun_ALL_Location.pdf"),width = 20,height = 12)


UPA <- dFrame[dFrame$Location %in% c("UPA","DSA") & dFrame$Class_Family %in% taxa,]
taxaPresent <- unique(UPA$Class_Family)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]
palette <- colors$hexID

UPA$Location <- factor(UPA$Location, levels = c("UPA","DSA"))

plot <- ggplot(UPA, aes(fill=Class_Family, y=Average, x=Location)) +
  geom_bar(stat="identity") +
  labs(x="", y="")+
  scale_fill_manual(values=palette)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(strip.text = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.text.y = element_text(size = 30))+
  scale_x_discrete(labels = c("UPA","DSA"))+
  scale_y_continuous(labels = scientific)

ggsave(plot, file=paste0(output, taxa_input, "_absAbun_ALL_Location_inset.pdf"),width = 20,height = 10)



