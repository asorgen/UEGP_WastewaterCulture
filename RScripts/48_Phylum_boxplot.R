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
inputModule = dir(pipeRoot, pattern="Gather", full.names=TRUE)
input = file.path(inputModule,"output/")
statModule = dir(pipeRoot, pattern="GLM_CULvENV", full.names=TRUE)
stat = file.path(statModule,"output/Phylum/")
output = file.path(moduleDir,"output/")

df=read.table(paste0(input, 'Phylum_table.tsv'),sep='\t',header = TRUE)
unique(df$Phylum)
df$Phylum = gsub("_","Other",df$Phylum)
df$Phylum=gsub("TM6..Dependentiae.","TM6 (Dependentiae)",df$Phylum)

glm <- read.table(paste0(stat, 'Phylum_GLM_Normalized_Combined.tsv'),sep='\t',header = TRUE)
glm <- glm[glm$pValuesSampleTypeAdjusted < 0.05, ]

plot <- ggplot(df, aes(y=Rel_Abun, x=Phylum, fill=SampleType)) +
  geom_boxplot()+
  labs(x="Phylum", y="Relative abundance")+
  # geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_minimal()+
  geom_signif(annotation=formatC("***", digits=1), y_position=1.05, xmin=1.75, xmax=2.25, tip_length = c(0, 0),textsize = 3)+
  geom_signif(annotation=formatC("***", digits=1), y_position=1.05, xmin=3.75, xmax=4.25, tip_length = c(0, 0),textsize = 3)+
  scale_fill_brewer(palette="Accent")+
  theme(axis.title.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18))+
  theme(axis.text.x = element_text(size = 10, angle = 75, hjust = 1))+
  theme(axis.text.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_text(size = 15))+
  # scale_x_discrete(labels = c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","No Antibiotic"))+
  theme(legend.position = "right")
ggsave(plot, file=paste(output, "Phylum_boxplot_relabun.pdf",sep = ""), width = 10, height = 7)
