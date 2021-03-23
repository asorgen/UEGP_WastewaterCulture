#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-05-21
#Description: Gathers normalized, relative abundance, absolute count data and merges into a single, large table.

library(gdata)
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)


rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
normModule = dir(pipeRoot, pattern="Normalize", full.names=TRUE)
normInput = file.path(normModule,"output/")
output = file.path(moduleDir,"output/")

input=c("Phylum","Class","Order","Family","Genus", "Species","Filtered_Family")
level = 1

for (j in 1:length(input)){
  inputname=paste0(normInput,input[j],"/",input[j])
  
  Normalized=read.table(file=paste0(inputname,"_Normalized.tsv"),sep="\t",header = TRUE)
  raw=read.table(file=paste0(inputname,"_raw.tsv"),sep="\t",header = TRUE)
  rel_abun=read.table(file=paste0(inputname,"_relabun.tsv"),sep="\t",header = TRUE)
  abs_abun=read.table(file=paste0(inputname,"_absabun.tsv"),sep="\t",header = TRUE)
  Normalized_CFU=read.table(file=paste0(inputname,"_Normalized_CFU.tsv"),sep="\t",header = TRUE)
  level = level + 1
  
  if (input[j] == "Filtered_Family") {
    level = 5
  }
  
  startTaxaIndex <- which(colnames(Normalized)=="CFU.mL")+1
  
  new <- gather(Normalized, "Taxon","Raw.Norm", startTaxaIndex:ncol(Normalized))
  new$Taxon=gsub(pattern ="D_[0-9]__",x = new$Taxon,replacement = "_")
  new$Taxon=gsub(pattern ="_Bacteria",x = new$Taxon,replacement = "Bacteria")
  other <- paste(replicate(level,"Other"), collapse = "._")
  new$Taxon=gsub(pattern ="Other",x = new$Taxon,replacement = other)
  string <- strsplit(as.character(new$Taxon),split = "._")
 
  temp_string=do.call(rbind,string)
  new <- cbind(temp_string,new)

  colnames(new)[colnames(new)=="1"] <- "Kingdom"
  colnames(new)[colnames(new)=="2"] <- "Phylum"
  colnames(new)[colnames(new)=="3"] <- "Class"
  colnames(new)[colnames(new)=="4"] <- "Order"
  colnames(new)[colnames(new)=="5"] <- "Family"
  colnames(new)[colnames(new)=="6"] <- "Genus"
  colnames(new)[colnames(new)=="7"] <- "Species"
  
  raw2 <- gather(raw, "Taxon","Raw_Counts", startTaxaIndex:ncol(raw))
  Raw_Counts <- raw2$Raw_Counts
  new <- cbind(new,Raw_Counts)
  
  rel_abun2 <- gather(rel_abun, "Taxon","Abundance", startTaxaIndex:ncol(rel_abun))
  Rel_Abun <- rel_abun2$Abundance
  new <- cbind(new,Rel_Abun)
  
  abs_abun2 <- gather(abs_abun, "Taxon","Abundance", startTaxaIndex:ncol(abs_abun))
  Abs_Abun <- abs_abun2$Abundance
  new <- cbind(new,Abs_Abun)
  
  Normalized_CFU2 <- gather(Normalized_CFU, "Taxon","Abundance", startTaxaIndex:ncol(Normalized_CFU))
  Norm.AbsAbun <- Normalized_CFU2$Abundance
  new <- cbind(new,Norm.AbsAbun)
  
  new <- new[ , -which(names(new) %in% c("Taxon"))]
  outputname=paste0(output,input[j])
  fwrite(x=new,file=paste0(outputname,"_table.tsv"),row.names = FALSE,sep="\t")
}


input=c("Filtered_Family")
level = 5

for (j in 1:length(input)){
  inputname=paste0(normInput,input[j],"/",input[j])
  
  Normalized=read.table(file=paste0(inputname,"_Normalized.tsv"),sep="\t",header = TRUE)
  raw=read.table(file=paste0(inputname,"_raw.tsv"),sep="\t",header = TRUE)
  rel_abun=read.table(file=paste0(inputname,"_relabun.tsv"),sep="\t",header = TRUE)
  abs_abun=read.table(file=paste0(inputname,"_absabun.tsv"),sep="\t",header = TRUE)
  Normalized_CFU=read.table(file=paste0(inputname,"_Normalized_CFU.tsv"),sep="\t",header = TRUE)
  
  new <- gather(Normalized, "Taxon","Raw.Norm", 17:ncol(Normalized))
  new$Taxon=gsub(pattern ="D_[0-9]__",x = new$Taxon,replacement = "_")
  new$Taxon=gsub(pattern ="_Bacteria",x = new$Taxon,replacement = "Bacteria")
  other <- paste(replicate(level,"Other"), collapse = "._")
  new$Taxon=gsub(pattern ="Other",x = new$Taxon,replacement = other)
  string <- strsplit(as.character(new$Taxon),split = "._")
  
  temp_string=do.call(rbind,string)
  new <- cbind(temp_string,new)
  
  colnames(new)[colnames(new)=="1"] <- "Kingdom"
  colnames(new)[colnames(new)=="2"] <- "Phylum"
  colnames(new)[colnames(new)=="3"] <- "Class"
  colnames(new)[colnames(new)=="4"] <- "Order"
  colnames(new)[colnames(new)=="5"] <- "Family"
  colnames(new)[colnames(new)=="6"] <- "Genus"
  
  raw2 <- gather(raw, "Taxon","Raw_Counts", 17:ncol(raw))
  Raw_Counts <- raw2$Raw_Counts
  new <- cbind(new,Raw_Counts)
  
  rel_abun2 <- gather(rel_abun, "Taxon","Abundance", 17:ncol(rel_abun))
  Rel_Abun <- rel_abun2$Abundance
  new <- cbind(new,Rel_Abun)
  
  abs_abun2 <- gather(abs_abun, "Taxon","Abundance", 17:ncol(abs_abun))
  Abs_Abun <- abs_abun2$Abundance
  new <- cbind(new,Abs_Abun)
  
  Normalized_CFU2 <- gather(Normalized_CFU, "Taxon","Abundance", 17:ncol(Normalized_CFU))
  Norm.AbsAbun <- Normalized_CFU2$Abundance
  new <- cbind(new,Norm.AbsAbun)
  
  new <- new[ , -which(names(new) %in% c("Taxon"))]
  outputname=paste0(output,input[j])
  fwrite(x=new,file=paste0(outputname,"_table.tsv"),row.names = FALSE,sep="\t")
}


