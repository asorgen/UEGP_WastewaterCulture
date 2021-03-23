#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description:

## Libraries
library(ggplot2)
library(grid)
library(scales)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Condition_RelAbun", full.names=TRUE)
input = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

taxa <- "Filtered_Family"

dFrame=read.table(paste0(input, taxa, "_Abundance_by_Condition.tsv"),sep='\t',header = TRUE)
dFrame$Class_Family = paste(dFrame$Class,dFrame$Family,sep = " - ")
dFrame2 <- dFrame[dFrame$Percentage >= 0.01,]
taxa <- unique(dFrame2$Class_Family)
dFrame <- dFrame[dFrame$Class_Family %in% taxa,]

dFrame$Class_Family=gsub("Other - Other","Other",dFrame$Class_Family)
dFrame$Class_Family=gsub("Bacilli - Family.XII","Bacilli - Family XII",dFrame$Class_Family)
dFrame$Class_Family=gsub("Flavobacteriia - NS9.marine.group","Flavobacteriia - NS9 marine group",dFrame$Class_Family)
dFrame$Class_Family=gsub("Sphingobacteriia - NS11.12.marine.group","Sphingobacteriia - NS11-12 marine group",dFrame$Class_Family)
dFrame$Class_Family=gsub("Cyanobacteria - FamilyI","Cyanobacteria - Family I",dFrame$Class_Family)

colorModule = dir(pipeRoot, pattern="PlotPrep", full.names=TRUE)
colorPath = file.path(colorModule,"output/")
colorIDs=read.table(paste0(colorPath,'colorIDs.tsv'),sep='\t',header = TRUE)
taxaPresent <- unique(dFrame$Class_Family)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]
palette <- colors$hexID

dFrame$Condition <- factor(dFrame$Condition, levels = c("Environ","BT_LB","RT_LB","BT_R2A","RT_R2A"))

Loc <- c("Environmental","LB agar at \n Body temp","LB agar at \n Room temp","R2A agar at \n Body temp","R2A agar at \n Room temp")
names(Loc) <- c("Environ","BT_LB","RT_LB","BT_R2A","RT_R2A")

plot2 <- ggplot(dFrame, aes(fill=Class_Family, y=Average, x=Condition)) +
  geom_bar(stat="identity",position = "fill") +
  labs(x="Condition", y="Relative Abundance")+
  scale_fill_manual(values=palette)+
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 22))+
  theme(strip.text = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.title.x = element_text(size = 30))+
  scale_x_discrete(labels = Loc)+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  guides(fill=guide_legend(title="Family",ncol = 1))+
  scale_y_continuous(labels = scientific)

ggsave(plot2, file=paste0(output, "Condition_relative_abundance.pdf"),width = 20,height = 15)
