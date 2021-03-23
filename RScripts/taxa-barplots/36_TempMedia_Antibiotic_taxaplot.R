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
inputModule = dir(pipeRoot, pattern="RelAbun_TempMedia_Antibiotic", full.names=TRUE)
input = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

taxa <- "Filtered_Family"

dFrame=read.table(paste0(input, taxa, "_Abundance_by_TempMedia_and_Antibiotic.tsv"),sep='\t',header = TRUE)

dFrame$Class_Family = paste(dFrame$Class,dFrame$Family,sep = " - ")
dFrame2 <- dFrame[dFrame$Percentage >= 0.01,]
taxa <- unique(dFrame2$Class_Family)
dFrame <- dFrame[dFrame$Class_Family %in% taxa,]

dFrame$Class_Family=gsub("Other - Other","Other",dFrame$Class_Family)
dFrame$Class_Family=gsub("Bacilli - Family.XII","Bacilli - Family XII",dFrame$Class_Family)

colorModule = dir(pipeRoot, pattern="PlotPrep", full.names=TRUE)
colorPath = file.path(colorModule,"output/")
colorIDs=read.table(paste0(colorPath,'colorIDs.tsv'),sep='\t',header = TRUE)
taxaPresent <- unique(dFrame$Class_Family)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]
palette <- colors$hexID

dFrame$Condition <- factor(dFrame$Condition, levels = c("BT_LB","RT_LB","BT_R2A","RT_R2A"))
dFrame$Antibiotic <- factor(dFrame$Antibiotic, levels = c("Amp","Cip","Dox","Sulf","Neg"))

ab <- c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","No Antibiotic")
names(ab) <- c("Amp","Cip","Dox","Sulf","Neg")
Loc <- c("LB agar at \n Body temp","LB agar at \n Room temp","R2A agar at \n Body temp","R2A agar at \n Room temp")
names(Loc) <- c("BT_LB","RT_LB","BT_R2A","RT_R2A")

plot <- ggplot(dFrame, aes(fill=Class_Family, y=Average, x=Antibiotic)) +
  geom_bar(stat="identity") +
  labs(x="Antibiotic", y="Average Abundance (CFU/mL)")+
  scale_fill_manual(values=palette)+
  theme_bw()+
  theme(legend.title = element_text(size = 25), legend.text = element_text(size = 15))+
  theme(strip.text = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  scale_x_discrete(labels = ab)+
  facet_grid(rows = vars(Condition), labeller = labeller(Condition=Loc))+
  guides(fill=guide_legend(title="Family",ncol = 1))+
  scale_y_continuous(labels = scientific)
ggsave(plot, file=paste0(output, "Condition_Antibiotic_absolute_abundance.pdf"),width = 20,height = 10)

plot2 <- ggplot(dFrame, aes(fill=Class_Family, y=Average, x=Antibiotic)) +
  geom_bar(stat="identity",position = "fill") +
  labs(x="Antibiotic", y="Relative Abundance")+
  scale_fill_manual(values=palette)+
  theme(legend.title = element_text(size = 25), legend.text = element_text(size = 18))+
  theme(strip.text = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  scale_x_discrete(labels = ab)+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  facet_grid(rows = vars(Condition), labeller = labeller(Condition=Loc))+
  guides(fill=guide_legend(title="Family",ncol = 1))+
  scale_y_continuous(labels = scientific)

ggsave(plot2, file=paste0(output, "Condition_Antibiotic_relative_abundance.pdf"),width = 20,height = 10)
