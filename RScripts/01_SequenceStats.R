#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: Provide a summary table for sequence data including sample size, diversity statistics, etc.

## Libraries
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(plyr)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="MetaUpdate", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

df=read.table(paste0(inputPath, "metaUpdate.tsv"),sep="\t",header = TRUE)
sampType <- unique(df$SampleType)

type <- vector()
stat <- vector()
value <- vector()
index <- 1

for (i in 1:length(sampType)) {
  df2 <- df[ df$SampleType %in% sampType[i] ,]

  type[index] <- sampType[i]
  stat[index] <- "SampleSize"
  value[index] <- nrow(df2)
  index <- index + 1
  
  for (j in c(13,14,15)) {
    
    type[index] <- sampType[i]
    stat[index] <- paste0(colnames(df2)[j], "_Average")
    value[index] <- mean(df2[,j])
    index <- index + 1
    
    type[index] <- sampType[i]
    stat[index] <- paste0(colnames(df2)[j], "_StDeviation")
    value[index] <- sd(df2[,j])
    index <- index + 1
    
    type[index] <- sampType[i]
    stat[index] <- paste0(colnames(df2)[j], "_Min")
    value[index] <- min(df2[,j])
    index <- index + 1
    
    type[index] <- sampType[i]
    stat[index] <- paste0(colnames(df2)[j], "_Max")
    value[index] <- max(df2[,j])
    index <- index + 1
    
    type[index] <- sampType[i]
    stat[index] <- paste0(colnames(df2)[j], "_Total")
    value[index] <- sum(df2[,j])
    index <- index + 1
    
  }
  
}
statSummary <- data.frame(type, stat, value)

lmStats <- lm(SequenceFrequency ~ SampleType, data = df)
anova <- anova(lmStats)
pValue <- anova$"Pr(>F)"[1]

statSummary <- rbind(statSummary, c("", "Seqlm", pValue))

lmStats <- lm(observed_otus ~ SampleType, data = df)
anova <- anova(lmStats)
pValue <- anova$"Pr(>F)"[1]

statSummary <- rbind(statSummary, c("", "OTUlm", pValue))

lmStats <- lm(shannon ~ SampleType, data = df)
anova <- anova(lmStats)
pValue <- anova$"Pr(>F)"[1]

statSummary <- rbind(statSummary, c("", "Shannonlm", pValue))

write.table(statSummary, file=paste0(output, "sequenceStatSummary.tsv"), sep="\t",row.names=FALSE)

