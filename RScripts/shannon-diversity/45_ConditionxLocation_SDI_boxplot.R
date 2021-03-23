#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: 

## Libraries
library(ggplot2)
library(grid)
library(scales)
library(ggsignif)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="MetaUpdate", full.names=TRUE)
input = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

df=read.table(paste0(input, "metaUpdate.tsv"),sep="\t",header = TRUE)

loc <- c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA")
df <- df[df$Location %in% loc,]
df$Location=factor(df$Location, levels = c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))
df$Condition=paste(df$Media,df$Temperature,sep = " agar at ")
condition_labs=c("LB agar at \n Body temp","LB agar at \n Room temp","R2A agar at \n Body temp","R2A agar at \n Room temp", "Environmental")
names(condition_labs) <- c("LB agar at BT","LB agar at RT","R2A agar at BT","R2A agar at RT","n/a agar at n/a")
palette <- c("#5599E7","#9B35E8","#DF6B7B","#9CF270","#EA9178","#C7F088","#E25047","#CBA7F2","#9DCBCB","#EDCDD0")


plot <- ggplot(df, aes(y=shannon, x=Location, fill=Location)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
  labs(x="Location", y="Shannon Diversity Index")+
  facet_grid(cols = vars(Condition),labeller = labeller(Condition=condition_labs))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_x_discrete(labels = loc)+
  scale_fill_manual(values = palette) +
  theme(strip.text = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 12))+
  theme(axis.text.y = element_text(size = 18))+
  theme(legend.text = element_text(size = 18))+
  theme(legend.title = element_text(size = 20))+
  theme(legend.position = "none")
ggsave(plot, file=paste(output, "ConditionxLocation_SDI_boxplot.pdf",sep = ""), width = 20, height = 7)


