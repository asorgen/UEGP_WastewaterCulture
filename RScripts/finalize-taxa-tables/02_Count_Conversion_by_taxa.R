#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-05-21
#Description: Combine all unassigned/ambiguous taxa and designate as 'Other'


library(gdata)
library(data.table)
library(stringr)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")

taxa_input=c("Phylum","Class","Order","Family","Genus","Species")

for (i in 2:7){
  # dir.create(paste0(output, taxa_input[i-1]), showWarnings = FALSE)
  inputname=paste0(input, "taxonomy-barplot/level-",i,".csv")
  myT=read.table(file=inputname,header = TRUE,sep = ",")
  names(myT)[names(myT)=="index"]="SampleID"
  endAbundanceIndex <- which(colnames(myT)=="SampleType")-1
  no_meta=myT[,1:endAbundanceIndex]
  
  if (taxa_input[i-1] == "Phylum") {
    list=colnames(no_meta)
    num=grep(x = list,pattern = "Unassigned.__")
    reduced=no_meta[,-num]
    fwrite(x=reduced,file=paste(output, taxa_input[i-1],"_No_Meta.tsv",sep=""),row.names = FALSE,sep="\t")
    
  } else {
    
    list=colnames(no_meta)
    list=list[list!="Unassigned.__"]
    num=grep(x = list,pattern = "\\.__")
    num=c(num, grep(x = list, pattern = "mbiguous"))
    num=c(num, grep(x = list, pattern = "unculture"))
    list_other=list[num]
    reduced=no_meta[,-num]
    reduced$Other=rowSums(no_meta[,num])
    
    fwrite(x=reduced,file=paste(output, taxa_input[i-1],"_No_Meta.tsv",sep=""),row.names = FALSE,sep="\t")
    
  }
}
  
sessionInfo()
