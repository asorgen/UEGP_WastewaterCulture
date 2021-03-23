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
inputModule = dir(pipeRoot, pattern="Antibiotic_RelAbun", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")
colorModule <- dir(pipeRoot, pattern="PlotPrep", full.names=TRUE)
colortPath = file.path(colorModule,"output/")



## .) Create abundance plot
taxa_input=c("Filtered_Family")
dFrame=read.table(paste0(inputPath, taxa_input,'_Abundance_by_Antibiotic.tsv'),sep='\t',header = TRUE)
dFrame$Class_Family = paste(dFrame$Class,dFrame$Family,sep = " - ")
dFrame$Class_Family=gsub("Other - Other","Other",dFrame$Class_Family)
dFrame$Class_Family=gsub("Bacilli - Family.XII","Bacilli - Family XII",dFrame$Class_Family)

# write.table(dFrame, file=paste0(output, taxa_input,"_Abundance_by_Antibiotic_Figs.tsv"),
#             sep = "\t", row.names = FALSE)

dFrame2 <- dFrame[dFrame$Percentage >= 0.01,]
taxa <- unique(dFrame2$Class_Family)
dFrame <- dFrame[dFrame$Class_Family %in% taxa,]



taxaPresent <- unique(dFrame$Class_Family)
colorIDs=read.table(paste0(colortPath,'colorIDs.tsv'),sep='\t',header = TRUE)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]

palette <- colors$hexID
dFrame$Antibiotic <- factor(dFrame$Antibiotic, levels = c("Amp", "Cip", "Dox", "Sulf", "Neg"))

plot <- ggplot(dFrame, aes(fill=Class_Family, y=Average, x=Antibiotic)) +
  geom_bar(stat="identity") +
  labs(x="Antibiotic", y="Average Abundance (CFU/mL)")+
  scale_fill_manual(values=palette)+
  theme_bw()+
  theme(legend.title = element_text(size = 25), legend.text = element_text(size = 22))+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 22))+
  theme(axis.text.y = element_text(size = 22))+
  scale_x_discrete(labels =c("Ampicillin", "Ciprofloxacin", "Doxycycline", "Sulfamethoxazole", "No Antibiotic"))+
  guides(fill=guide_legend(title="Class - Family",ncol = 1))+
  scale_y_continuous(labels = scientific)

ggsave(plot, file=paste0(output, taxa_input, "_absAbun_Antibiotic.pdf"),width = 20,height = 10)




## .) Ciprofloxacin inset
Cip_inset <- dFrame[dFrame$Antibiotic %in% c("Cip") & dFrame$Class_Family %in% taxa,]
taxaPresent <- unique(Cip_inset$Class_Family)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]
palette <- colors$hexID

Cip_inset$Antibiotic <- factor(Cip_inset$Antibiotic)

plot <- ggplot(Cip_inset, aes(fill=Class_Family, y=Average, x=Antibiotic)) +
  geom_bar(stat="identity") +
  labs(x="", y="")+
  scale_fill_manual(values=palette)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(strip.text = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.text.y = element_text(size = 30))+
  scale_x_discrete(labels = c("Ciprofloxacin"))+
  scale_y_continuous(labels = scientific)

ggsave(plot, file=paste0(output, taxa_input, "_absAbun_Cip_inset.pdf"),width = 5,height = 10)




## .) Doxycycline inset
Dox_inset <- dFrame[dFrame$Antibiotic %in% c("Dox") & dFrame$Class_Family %in% taxa,]
taxaPresent <- unique(Dox_inset$Class_Family)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]
palette <- colors$hexID

Dox_inset$Antibiotic <- factor(Dox_inset$Antibiotic)

plot <- ggplot(Dox_inset, aes(fill=Class_Family, y=Average, x=Antibiotic)) +
  geom_bar(stat="identity") +
  labs(x="", y="")+
  scale_fill_manual(values=palette)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(strip.text = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.text.y = element_text(size = 30))+
  scale_x_discrete(labels = c("Doxycycline"))+
  scale_y_continuous(labels = scientific)

ggsave(plot, file=paste0(output, taxa_input, "_absAbun_Dox_inset.pdf"),width = 5,height = 10)





## .) Sulfamethoxazole inset
Sulf_inset <- dFrame[dFrame$Antibiotic %in% c("Sulf") & dFrame$Class_Family %in% taxa,]
taxaPresent <- unique(Sulf_inset$Class_Family)
colors <- colorIDs[colorIDs$taxaNames %in% taxaPresent,]
palette <- colors$hexID

Sulf_inset$Antibiotic <- factor(Sulf_inset$Antibiotic)

plot <- ggplot(Sulf_inset, aes(fill=Class_Family, y=Average, x=Antibiotic)) +
  geom_bar(stat="identity") +
  labs(x="", y="")+
  scale_fill_manual(values=palette)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(strip.text = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 30))+
  scale_x_discrete(labels = c("Sulfamethoxazole"))+
  scale_y_continuous(labels = scientific)

ggsave(plot, file=paste0(output, taxa_input, "_absAbun_Sulf_inset.pdf"),width = 5,height = 10)


