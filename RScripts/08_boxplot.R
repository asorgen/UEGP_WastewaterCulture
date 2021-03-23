#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-05-21
#Description: Creating average environmental OTU abundance found in cultures

rm(list=ls())

## 1.) Set libraries
message("1.) Set libraries")
#library('ggsignif')
library('ggplot2')
library('stringr')
message("1.) Set libraries - DONE!")



pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="OTUtable", full.names=TRUE)
input = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")


## 4.) Set up for boxplot
message("4.) Set up for boxplot")
df <- read.table(paste0(input, "percentage_table.tsv"), header = TRUE, sep = "\t")
df$x.axis <- "Culture"
df[is.na(df)] <- 0
meanVal <- aggregate(Env_Percentage ~  x.axis, df, mean)
message("4.) Set up for boxplot - DONE!")



## 5.) Generate boxplot
message("5.) Generate boxplot")
plot <- ggplot(df, aes(y=Env_Percentage, x=x.axis)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+
  stat_summary(fun = "mean", color="blue", geom="point", shape=18, position=position_dodge(.73), size=5)+
  labs(x="Cultured\nEnvironmental OTUs", y="Relative Abundance")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_bw()+
  scale_x_discrete(labels = c(" "))+
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 13))+
  geom_hline(yintercept=0.01, linetype="dashed", color = "red")+
  theme(legend.position = "right")

message("5.) Generate boxplot - DONE!")



## 6.) Save boxplot
message("6.) Save boxplot")
ggsave(plot, file=paste0(output, "OTU_boxplot.pdf"), width = 2.5, height = 7)
message("6.) Save boxplot - DONE!")

