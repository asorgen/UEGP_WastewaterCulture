rm(list=ls())
# library(plyr)
library(tidyverse)
library(ggplot2)
# library(pdp)

########## Making combined OTU table ##########
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="OTUtable", full.names=TRUE)
input = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")
inputPath = file.path(pipeRoot,"input/")

df <- read.table(paste0(input, "relabun_OTUtable.tsv"),sep = "\t",header = TRUE)
myT <- read.table(paste0(input, "raw_OTUtable.tsv"),sep = "\t",header = TRUE)
Edf <- df[df$SampleType=="Environ",] # create environmental sample only table
Cdf <- df[df$SampleType=="Culture",] # create culture sample only table
EmyT <- myT[myT$SampleType=="Environ",] # create environmental sample only table
CmyT <- myT[myT$SampleType=="Culture",] # create culture sample only table

CulAvgRelAbun <- vector()
EnvAvgRelAbun <- vector()
OTU <- vector()
CulOTUFreq <- vector()
EnvOTUFreq <- vector()
tTest.pVal <- vector()
CulSampleNum <- vector()
EnvSampleNum <- vector()
index <- 1
OTUstart <- which(colnames(df)=="SequenceFrequency")+1

for (i in OTUstart:ncol(df)) {
  OTU[index] <- colnames(df)[i]
  
  Cdf2 <- Cdf[Cdf[,i] > 0,]
  if (sum(Cdf2[,i]) > 0) {
    CulAvgRelAbun[index] <- mean(Cdf2[,i])
    CulSampleNum[index] <- nrow(Cdf2)
    CulOTUFreq[index] <- sum(CmyT[,i])
  } else {
    CulAvgRelAbun[index] <- 0
    CulSampleNum[index] <- 0
    CulOTUFreq[index] <- 0
  }
  
  Edf2 <- Edf[Edf[,i] > 0,]
  if (sum(Edf2[,i]) > 0) {
    EnvAvgRelAbun[index] <- mean(Edf2[,i])
    EnvSampleNum[index] <- nrow(Edf2)
    EnvOTUFreq[index] <- sum(EmyT[,i])
  } else {
    EnvAvgRelAbun[index] <- 0
    EnvSampleNum[index] <- 0
    EnvOTUFreq[index] <- 0
  }
  
  df2 <- df[df[,i] > 0,]
  if ((nrow(Cdf2) > 1) && (nrow(Edf2) > 1)) {
    t_test <- t.test(df2[,i] ~ df2$SampleType)
    tTest.pVal[index] <- t_test$p.value
  } else {
    tTest.pVal[index] <- 2
  }
  

  index <- index + 1
}

frame <- cbind(OTU, CulAvgRelAbun, EnvAvgRelAbun, CulSampleNum, EnvSampleNum, CulOTUFreq, EnvOTUFreq, tTest.pVal)
write.table(frame, file=paste0(output, "avg_OTU_relabun.tsv"), sep="\t",row.names=FALSE)

myT <- read.table(paste0(output, "avg_OTU_relabun.tsv"),sep = "\t",header = TRUE)
myT$OTU = gsub("X","",myT$OTU)

tax <- read.table(paste0(inputPath, "OTU-taxonomy2.txt"), sep = "\t", header = TRUE)

merge <- merge(tax, myT, by = "OTU")
write.table(merge,file=paste0(output, "merge.tsv"),sep="\t",row.names = FALSE)

myT2 <- merge[merge$CulSampleNum > 0,]
myT2 <- myT2[myT2$EnvSampleNum > 0,]

cols <- c("red", "orange",  "yellow", "green", "blue", "purple","violet", "pink", "black")
shapes <- c(16, 16, 16, 16, 16, 16, 1, 16, 16)

p <- which(myT2$Phylum %in% "Proteobacteria")
o <- which(!(myT2$Phylum %in% "Proteobacteria"))

plot <- ggplot()+
  geom_point(data = myT2, aes(x = CulAvgRelAbun, y = EnvAvgRelAbun, color = Phylum), size = 3)+
  theme_bw()+
  labs(x="Cultured OTU Relative Abundance", y="Environmental OTU Relative Abundance")+
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 18))+
  theme(axis.text.y = element_text(size = 18))+
  # scale_color_manual(values = cols)
  xlim(0,1)+
  ylim(0,1)+
  theme(legend.title = element_text(size = 18), legend.text = element_text(size = 15))+
  geom_abline(intercept = 0, slope = 1, col = "red")
ggsave(plot, file=paste0(output, "Phylum_OTU_scatterplot.pdf"), width = 12, height = 10)

plot <- ggplot()+
  geom_point(data = myT2, aes(x = CulAvgRelAbun, y = EnvAvgRelAbun, color = Class), size = 2)+
  theme_bw()+
  labs(x="Culture OTU Relative Abundance", y="Environmental OTU Relative Abundance")+
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.title.x = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 18))+
  theme(axis.text.y = element_text(size = 18))+
  xlim(0,1)+
  ylim(0,1)+
  geom_abline(intercept = 0, slope = 1, col = "red")


ggsave(plot, file=paste0(output, "Class_OTU_scatterplot.pdf"), width = 10, height = 10)
