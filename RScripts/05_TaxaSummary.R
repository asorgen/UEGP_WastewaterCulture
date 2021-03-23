#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-25-21
#Description: 

## Libraries
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(plyr)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Gather", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

column <- 1
taxa=c("Phylum","Class","Order","Family","Genus")
SampleType <- vector()
TaxaLevel <- vector()
Count <- vector()
index <- 1

for (i in 1:length(taxa)) {
    df <- read.table(paste0(inputPath, taxa[i],'_table.tsv'),sep='\t',header = TRUE)
    df <- df[df$Raw_Counts > 0,]
    column <- column + 1
    df <- df[!(df[,column] %in% c("_", "Other")),]
    
  
    for (j in c("Culture","Environ")) {
      SampleType[index] <- j
      df2 <- df[df$SampleType %in% j, ]
      TaxaLevel[index] <- colnames(df2)[column]
      Count[index] <- length(unique(df2[,column]))
      index <- index +1
    }
}

dFrame <- data.frame(SampleType, TaxaLevel, Count)

write.table(dFrame, file=paste0(output, "TaxaSummary.tsv"), sep="\t",row.names=FALSE)
