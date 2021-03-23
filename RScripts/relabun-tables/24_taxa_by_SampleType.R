#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-07-21
#Description: 

rm(list=ls())
taxa=c("Phylum","Class","Order","Family","Genus","Filtered_Family")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Gather", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

L <- 1

for (z in 1:length(taxa)) {
  frame=read.table(paste0(inputPath, taxa[z],'_table.tsv'),sep='\t',header = TRUE)
  L <- L + 1
  frame$Taxa <- paste(frame$Kingdom,frame$Phylum,frame$Class,frame$Order,frame$Family,frame$Genus,sep = ">")
  list <- c("Culture","Environ")
  
  if (taxa[z] == "Filtered_Family") {
    L <- 5
  }
  
  ## SampleType averages
  dFrame=data.frame()
  for (x in 1:length(list)) {
    SampleType=paste(list[x])
    df=frame[frame$SampleType %in% SampleType,]
    
    TAXA=factor(unique(frame$Taxa))
    Total_CFU=sum(df$Raw_Counts)
    
    for (i in 1:length(TAXA)) {
      Taxa=paste(TAXA[i])
      df3=df[df$Taxa %in% Taxa,]
      Taxa_Total=sum(df3$Raw_Counts)
      Average=mean(df3$Raw_Counts)
      SD=sd(df3$Raw_Counts)
      Percentage=Taxa_Total / Total_CFU
      addition=data.frame(SampleType,Taxa,Taxa_Total,Total_CFU,Percentage,Average,SD)
      dFrame=rbind(dFrame,addition)
    }
  }
  string <- strsplit(as.character(dFrame$Taxa),split = ">")
  temp_string=do.call(rbind,string)
  dFrame <- cbind(temp_string[,1:L],dFrame)
  dFrame = subset(dFrame, select = -c(Taxa) )
  
  colnames(dFrame)[colnames(dFrame)=="1"] <- "Kingdom"
  colnames(dFrame)[colnames(dFrame)=="2"] <- "Phylum"
  colnames(dFrame)[colnames(dFrame)=="3"] <- "Class"
  colnames(dFrame)[colnames(dFrame)=="4"] <- "Order"
  colnames(dFrame)[colnames(dFrame)=="5"] <- "Family"
  colnames(dFrame)[colnames(dFrame)=="6"] <- "Genus"
  colnames(dFrame)[colnames(dFrame)=="7"] <- "Species"
  
  write.table(dFrame, file=paste0(output, taxa[z],"_Abundance_by_SampleType.tsv"),sep = "\t", row.names = FALSE)
  
}

