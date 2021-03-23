#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-05-21
#Description: 

rm(list=ls())
# library(plyr)
library(tidyverse)

########## Making combined OTU table ##########
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input <- paste0(pipeRoot, "/input/")
output = file.path(moduleDir,"output/")

c1 <- read.table(paste0(input, "OR-OTUfeature-table.tsv"),sep = "\t",header = TRUE)
c2 <- t(c1)
colnames(c2) <- c2[1,]
c2 <- c2[-1,]

SampleID <- rownames(c2)
c3 <- cbind(SampleID, c2)

metaModule = dir(pipeRoot, pattern="MetaUpdate", full.names=TRUE)
metaInput = file.path(metaModule,"output/")
meta <- read.table(paste0(metaInput, "metaUpdate.tsv"),sep = "\t",header = TRUE)
reduced=merge(meta,c3,by="SampleID")
write.table(reduced, file=paste0(output, "raw_OTUtable.tsv"), sep="\t",row.names=FALSE)


reduced <- read.table(paste0(output, "raw_OTUtable.tsv"), sep="\t",header = TRUE)
rel_abun=reduced
OTUstart <- which(colnames(reduced)=="SequenceFrequency")+1
for (i in 1:nrow(reduced)){
  rel_abun[i,OTUstart:ncol(reduced)] = reduced[i,OTUstart:ncol(reduced)] / rowSums(reduced[i,OTUstart:ncol(reduced)])
}

write.table(rel_abun, file=paste0(output, "relabun_OTUtable.tsv"), sep="\t",row.names=FALSE)


Normalized=reduced
for (i in 1:nrow(reduced)){
  Normalized[i,OTUstart:ncol(reduced)] = 
    log10(((reduced[i,OTUstart:ncol(reduced)]) / (rowSums(reduced[,OTUstart:ncol(reduced)])[i]) * 
             rowMeans(reduced[,OTUstart:ncol(reduced)])[i]) + 1)
}
write.table(x=Normalized,file=paste0(output,"Normalized_OTUtable.tsv"),row.names = FALSE,sep="\t")


new <- gather(reduced, "OTU","Count", OTUstart:ncol(reduced))
df <- gather(rel_abun, "OTU","Rel.Abun", OTUstart:ncol(rel_abun))
Rel.Abun <- df$Rel.Abun
new <- cbind(new,Rel.Abun)

write.table(new, file=paste0(output, "OTUtable.tsv"), sep="\t",row.names=FALSE)


########## ENVIRONMENTAL PERCENTAGE TABLE: ##########

df <- read.table(paste0(output, "OTUtable.tsv"),sep = "\t",header = TRUE)

dFrame <- data.frame()
S <- c("Mallard", "Sugar")

for (s in 1:length(S)) {
  Site <- paste(S[s])
  myS <- df[df$Site %in% Site,]
  L <- unique(myS$Location)
  
  for (l in 1:length(L)) {
    Location <- paste(L[l])
    myL <- myS[myS$Location %in% Location,]
    myL <- myL[myL$Count > 0,]
    myC <- myL[myL$SampleType=="Culture",]
    sampleOTUs <- unique(myC$OTU)
    numCulOTUs <- length(unique(myC$OTU))
    
    ENV <- df[df$SampleType=="Environ",]
    ENV <- ENV[ENV$Site %in% Site & ENV$Location %in% Location,]
    ENV <- ENV[ENV$Count > 0,]
    sampleID <- unique(ENV$SampleID)
    
    for (x in 1:length(sampleID)) {
      SampleID <- paste(sampleID[x])
      df4 <- ENV[ENV$SampleID %in% SampleID,]
      eSampleTotal <- sum(df4$Count)
      numEnvOTUs <- length(unique(df4$OTU))
      ENV2 <- df4[df4$OTU %in% sampleOTUs,]
      num_of_culturedEnvOTUs <- length(unique(ENV2$OTU))
      culturedEnvTotal <- sum(ENV2$Count)
      Env_Percentage <- culturedEnvTotal/eSampleTotal
      frame <- cbind(Site,	Location,	SampleID,	numCulOTUs,	eSampleTotal,	numEnvOTUs,	num_of_culturedEnvOTUs,	culturedEnvTotal,	Env_Percentage)
      dFrame <- rbind(dFrame,frame)
    }
  }
}

write.table(dFrame, file=paste0(output, "percentage_table.tsv"), sep="\t",row.names=FALSE)

