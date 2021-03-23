#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-05-21
#Description: Normalize count data for each taxonomic OTU table


library(gdata)
library(data.table)
library(stringr)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
metaModule = dir(pipeRoot, pattern="MetaUpdate", full.names=TRUE)
metaInput = file.path(metaModule,"output/")
countModule = dir(pipeRoot, pattern="CountConvert", full.names=TRUE)
countInput = file.path(countModule,"output/")
output = file.path(moduleDir,"output/")

taxa_input=c("Phylum","Class","Order","Family","Filtered_Family","Genus","Species")


##### Filter out taxonomic families present in less than 3% relative abundance
meta_data=read.table(paste0(metaInput, "metaUpdate.tsv"),sep="\t",header=TRUE)

df = read.table(paste0(countInput, "Family_No_Meta.tsv"), header = TRUE, sep = "\t")
duplicates=c("C223","C179","C143","C119")
meta_data=meta_data[!(meta_data$SampleID %in% duplicates),]
meta_data$CFU.mL = ((meta_data$Counts/50)*1000)*meta_data$Dilution
reduced=merge(meta_data,df,by="SampleID")
startTaxaIndex <- which(colnames(reduced)=="CFU.mL")+1

rel_abun=reduced
for (x in 1:nrow(reduced)){
  rel_abun[x,startTaxaIndex:ncol(reduced)]=((reduced[x,startTaxaIndex:ncol(reduced)])/(rowSums(reduced[,startTaxaIndex:ncol(reduced)])[x]))
}

colMax <- function(data) sapply(data, max, na.rm = TRUE)
list = colMax(rel_abun[,startTaxaIndex:ncol(rel_abun)])
list=list[list < 0.03]
rare = names(list)

new=reduced[,!colnames(reduced) %in% rare]
out=reduced[,colnames(reduced) %in% rare]
new$Other = rowSums(out) + new$Other
startMeta <- which(colnames(reduced)=="SampleType")
endMeta <- which(colnames(reduced)=="CFU.mL")
new=new[,-(startMeta:endMeta)]

fwrite(x=new,file=paste0(countInput, "Filtered_Family_No_Meta.tsv"),row.names = FALSE,sep="\t")


for (j in 1:length(taxa_input)) {
  dir.create(paste0(output, taxa_input[j]), showWarnings = FALSE)
  taxaDir <- paste0(output, taxa_input[j], "/")
  inputname=paste0(countInput, taxa_input[j])
  myT=read.table(file=paste0(inputname, "_No_Meta.tsv"),header = TRUE,sep = "\t")
  meta_data=read.table(paste0(metaInput, "metaUpdate.tsv"),sep="\t",header=TRUE)
  duplicates=c("C223","C179","C143","C119")
  meta_data=meta_data[!(meta_data$SampleID %in% duplicates),]
  meta_data$CFU.mL = ((meta_data$Counts/50)*1000)*meta_data$Dilution
  
  reduced=merge(meta_data,myT,by="SampleID")
  outputname=paste0(taxaDir, taxa_input[j])
  fwrite(x=reduced,file=paste0(outputname,"_raw.tsv"),row.names = FALSE,sep="\t")
  
  startTaxaIndex <- which(colnames(reduced)=="CFU.mL")+1
  rel_abun=reduced
  for (x in 1:nrow(reduced)){
    rel_abun[x,startTaxaIndex:ncol(reduced)]=((reduced[x,startTaxaIndex:ncol(reduced)])/(rowSums(reduced[,startTaxaIndex:ncol(reduced)])[x]))
  }
  fwrite(x=rel_abun,file=paste0(outputname,"_relabun.tsv"),row.names = FALSE,sep="\t")
  
  abs_abun=rel_abun
  for (x in 1:nrow(rel_abun)){
    abs_abun[x,startTaxaIndex:ncol(rel_abun)] = rel_abun[x,startTaxaIndex:ncol(rel_abun)] * rel_abun$CFU.mL[x]
  }
  fwrite(x=abs_abun,file=paste0(outputname,"_absabun.tsv"),row.names = FALSE,sep="\t")
  
  log_abun=abs_abun
  for (x in 1:nrow(abs_abun)){
    log_abun[x,startTaxaIndex:ncol(abs_abun)] = log10(abs_abun[x,startTaxaIndex:ncol(abs_abun)]+1) 
  }
  fwrite(x=log_abun,file=paste0(outputname,"_log10_absabun.tsv"),row.names = FALSE,sep="\t")
  
  
  Normalized=reduced
  for (i in 1:nrow(reduced)){
    Normalized[i,startTaxaIndex:ncol(reduced)] = 
      log10(((reduced[i,startTaxaIndex:ncol(reduced)]) / (rowSums(reduced[,startTaxaIndex:ncol(reduced)])[i]) * 
               rowMeans(reduced[,startTaxaIndex:ncol(reduced)])[i]) + 1)
  }
  fwrite(x=Normalized,file=paste0(outputname,"_Normalized.tsv"),row.names = FALSE,sep="\t")
  
  
  Norm_CFU=abs_abun
  for (i in 1:nrow(abs_abun)){
    Norm_CFU[i,startTaxaIndex:ncol(abs_abun)] = 
      log10(((abs_abun[i,startTaxaIndex:ncol(abs_abun)]) / (rowSums(abs_abun[,startTaxaIndex:ncol(abs_abun)])[i]) * 
               rowMeans(abs_abun[,startTaxaIndex:ncol(abs_abun)])[i]) + 1)
  }
  fwrite(x=Norm_CFU,file=paste0(outputname,"_Normalized_CFU.tsv"),row.names = FALSE,sep="\t")
  
}
