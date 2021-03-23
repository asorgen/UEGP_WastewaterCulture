#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-05-21
#Description: 

rm(list=ls())

## Libraries
# library(data.table)
# library(stringr)
library(tidyr)
library(ggplot2)

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Gather", full.names=TRUE)
input = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")





## 1.) Find average normalized counts for each OTU per sample site
taxa_input=c("Phylum","Class","Order","Family","Genus","Filtered_Family")

L <- 1

for (z in 1:length(taxa_input)) {
  frame=read.table(paste(input, taxa_input[z],'_table.tsv',sep=""),sep='\t',header = TRUE)
  frame <- frame[frame$SampleType == "Culture",]
  L <- L + 1
  site <- c("Mallard","Sugar")

  if (taxa_input[z] == "Filtered_Family") {
    L <- 5
  }
  
  ## Site averages
  dFrame=data.frame()
  for (x in 1:length(site)) {
    Site=paste(site[x])
    df=frame[frame$Site %in% Site,]
    
    TAXA=factor(unique(frame[,L]))
    
    for (i in 1:length(TAXA)) {
      Taxa=paste(TAXA[i])
      df3=df[df[,L] %in% Taxa,]
      Average=mean(df3$Norm.AbsAbun)
      addition=data.frame(Site,Taxa,Average)
      dFrame=rbind(dFrame,addition)
    }
  }
  write.table(dFrame, file=paste(output,taxa_input[z],"_by_Site_ALL_abs_abun_normalized.tsv",sep = ""),sep = "\t", row.names = FALSE)
}





## 2.) Find average normalized counts for each OTU per sample site for ARB only
L <- 1

for (z in 1:length(taxa_input)) {
  frame=read.table(paste(input, taxa_input[z],'_table.tsv',sep=""),sep='\t',header = TRUE)
  frame <- frame[frame$SampleType == "Culture",]
  frame <- frame[!(frame$Antibiotic == "Neg"),]
  L <- L + 1
  site <- c("Mallard","Sugar")
  
  if (taxa_input[z] == "Filtered_Family") {
    L <- 5
  }
  
  ## Site averages
  dFrame=data.frame()
  for (x in 1:length(site)) {
    Site=paste(site[x])
    df=frame[frame$Site %in% Site,]
    
    TAXA=factor(unique(frame[,L]))
    
    for (i in 1:length(TAXA)) {
      Taxa=paste(TAXA[i])
      df3=df[df[,L] %in% Taxa,]
      Average=mean(df3$Norm.AbsAbun)
      addition=data.frame(Site,Taxa,Average)
      dFrame=rbind(dFrame,addition)
    }
  }
  write.table(dFrame, file=paste(output,taxa_input[z],"_by_Site_ARB_abs_abun_normalized.tsv",sep = ""),sep = "\t", row.names = FALSE)
}







## 3.) Generate scatter plot for abundance of OTUs per sampling site
statinputModule = dir(pipeRoot, pattern="GLM_CFU_ALL", full.names=TRUE)
statinput = file.path(statinputModule,"output/")
taxa_input=c("Filtered_Family")

dFrame=read.table(paste(output,taxa_input,'_by_Site_ALL_abs_abun_normalized.tsv',sep=""),sep='\t',header = TRUE)
stats <- read.table(paste(statinput,taxa_input,"/",taxa_input,"_GLM_NormalizedCFU_Culture.tsv",sep = ""), sep = '\t', header = TRUE)
stats$Name=gsub(pattern ="D_[0-9]__",x = stats$Name,replacement = "_")
stats$Name=gsub(pattern ="_Bacteria",x = stats$Name,replacement = "Bacteria")
string <- strsplit(as.character(stats$Name),split = "._")
temp_string=do.call(rbind,string)
stats <- cbind(temp_string,stats)
there <- paste(stats[,5])
sig <- stats[stats$pValuesSiteAdjusted < 0.05,]
sig <- sig[!(sig$Name=="Other"),]
marks <- paste(sig[,5])

data_spread <- dFrame %>%
  spread(key = Site, value = Average)
h <- which(data_spread$Taxa %in% marks)
n <- which(!(data_spread$Taxa %in% marks) & !(data_spread$Taxa %in% there))

plot <- ggplot()+
  geom_point(data = data_spread, aes(Mallard,Sugar), color = "black", size = 2)+
  geom_point(data = data_spread[h, ], aes(Mallard,Sugar), color = "red", size = 2)+
  geom_point(data = data_spread[n, ], aes(Mallard,Sugar), color = "grey", size = 2)+
  theme(legend.position = "none")+
  geom_text(data= data_spread[h, ],aes(Mallard,Sugar,label=Taxa), color = "red", hjust=-0.07, size=3)+
  # geom_text(label=data_spread$Taxa, nudge_x = .3,nudge_y = .3, check_overlap = TRUE)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 17))+
  theme(axis.title.y = element_text(size = 17))+
  geom_abline(intercept = 0, slope = 1)

ggsave(plot, file=paste(output, taxa_input,"_by_Site_ALL_scatter.pdf",sep = ""), width = 7, height = 7)







## 4.) Generate scatter plot for abundance of OTUs per sampling site for ARB only
statinputModule = dir(pipeRoot, pattern="GLM_CFU_ARB", full.names=TRUE)
statinput = file.path(statinputModule,"output/")

dFrame=read.table(paste(output,taxa_input,'_by_Site_ARB_abs_abun_normalized.tsv',sep=""),sep='\t',header = TRUE)
stats <- read.table(paste(statinput,taxa_input,"/",taxa_input,"_GLM_NormalizedCFU_Culture_ARB_only.tsv",sep = ""), sep = '\t', header = TRUE)
stats$Name=gsub(pattern ="D_[0-9]__",x = stats$Name,replacement = "_")
stats$Name=gsub(pattern ="_Bacteria",x = stats$Name,replacement = "Bacteria")
string <- strsplit(as.character(stats$Name),split = "._")
temp_string=do.call(rbind,string)
stats <- cbind(temp_string,stats)
there <- paste(stats[,5])
sig <- stats[stats$pValuesSiteAdjusted < 0.05,]
sig <- sig[!(sig$Name=="Other"),]
marks <- paste(sig[,5])

data_spread <- dFrame %>%
  spread(key = Site, value = Average)
h <- which(data_spread$Taxa %in% marks)
n <- which(!(data_spread$Taxa %in% marks) & !(data_spread$Taxa %in% there))

plot <- ggplot()+
  geom_point(data = data_spread, aes(Mallard,Sugar), color = "black", size = 2)+
  geom_point(data = data_spread[h, ], aes(Mallard,Sugar), color = "red", size = 2)+
  geom_point(data = data_spread[n, ], aes(Mallard,Sugar), color = "grey", size = 2)+
  theme(legend.position = "none")+
  geom_text(data= data_spread[h, ],aes(Mallard,Sugar,label=Taxa), color = "red", hjust=-0.07, size=3)+
  # geom_text(label=data_spread$Taxa, nudge_x = .3,nudge_y = .3, check_overlap = TRUE)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 17))+
  theme(axis.title.y = element_text(size = 17))+
  geom_abline(intercept = 0, slope = 1)

ggsave(plot, file=paste(output, taxa_input,"_by_Site_ARB_scatter.pdf",sep = ""), width = 7, height = 7)
