#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-25-21
#Description: Linear model comparisions between individual locations.

library(vegan)

rm(list=ls())
input=c("Filtered_Family","Phylum","Class","Order","Family","Genus")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Normalize", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")
args <- commandArgs(trailingOnly = TRUE)

loc1 <- args[1]
loc2 <- args[2]

for (z in 1:length(input)) {
  locations <- c(loc1, loc2)
  inputDir <- paste0(inputPath, input[z], "/")
  myT=read.table(paste0(inputDir, input[z],'_Normalized_CFU.tsv'),sep='\t',header = TRUE)
  myT <- myT[myT$SampleType == "Culture",]
  myT <- myT[myT$Location %in% locations,]
  myT <- myT[!(myT$Antibiotic == "Neg"),]
  
  pValuesLocation <- vector()
  Name <- vector()
  index <- 1
  indexes <- vector()
  absolute_indexes=vector()
  taxaStart <- which(colnames(myT)=="CFU.mL")+1
  
  
  for(i in taxaStart:ncol(myT)){
    if(sum(myT[,i]!=0) > nrow(myT)*0.1 ) 
    {
      Location <-factor(myT$Location)

      myLm <- lm(myT[,i] ~ Location)
      myAnova <- anova(myLm)
      pValuesLocation[index] <- myAnova$"Pr(>F)"[1]

      Name[index] <- names(myT)[i]
      indexes[index] <- index
      absolute_indexes[index]=i
      
      index <- index + 1
    }
  }
  # Making a dataframe out of all the p values and indexes
  dFrame <- data.frame(Name, indexes, absolute_indexes,pValuesLocation)
  # Performing BH adjustments on p values
  dFrame$pValuesLocationAdjusted<- p.adjust( dFrame$pValuesLocation, method = "BH")

  dFrame$Name=gsub(pattern ="D_[0-9]__",x = dFrame$Name,replacement = "_")
  dFrame$Name=gsub(pattern ="_Bacteria",x = dFrame$Name,replacement = "Bacteria")
  string <- strsplit(as.character(dFrame$Name),split = "._")
  
  temp_string=do.call(rbind,string)
  dFrame <- cbind(temp_string,dFrame)
  
  colnames(dFrame)[colnames(dFrame)=="1"] <- "Kingdom"
  colnames(dFrame)[colnames(dFrame)=="2"] <- "Phylum"
  colnames(dFrame)[colnames(dFrame)=="3"] <- "Class"
  colnames(dFrame)[colnames(dFrame)=="4"] <- "Order"
  colnames(dFrame)[colnames(dFrame)=="5"] <- "Family"
  colnames(dFrame)[colnames(dFrame)=="6"] <- "Genus"
  colnames(dFrame)[colnames(dFrame)=="7"] <- "Species"
  
  dir.create(paste0(output, input[z]), showWarnings = FALSE)
  outputDir <- paste0(output, input[z],"/")
  dir.create(paste0(outputDir,"PDFs/"), showWarnings = FALSE)
  PDFdir <- paste0(outputDir, "PDFs/")
  
  # Saving results as a text file
  write.table(dFrame, file=paste0(outputDir, input[z],"_GLM_NormalizedCFU_Culture_",loc1,"v",loc2,"_ARBonly.tsv"), sep="\t",row.names=FALSE)
  
  if (min(dFrame$pValuesLocationAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesLocationAdjusted < 0.05,]
    Location <- paste(loc1, loc2, sep = " - ")
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesLocationAdjusted
    dFrame2 <- cbind(Taxa, Location, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ARB.tsv"), sep="\t",row.names=FALSE)
  }
  
}


for (z in 1:length(input)) {
  locations <- c(loc1, loc2)
  inputDir <- paste0(inputPath, input[z], "/")
  myT=read.table(paste0(inputDir, input[z],'_Normalized_CFU.tsv'),sep='\t',header = TRUE)
  myT <- myT[myT$SampleType == "Culture",]
  myT <- myT[myT$Location %in% locations,]

  pValuesLocation <- vector()
  Name <- vector()
  index <- 1
  indexes <- vector()
  absolute_indexes=vector()
  taxaStart <- which(colnames(myT)=="CFU.mL")+1
  
  
  for(i in taxaStart:ncol(myT)){
    if(sum(myT[,i]!=0) > nrow(myT)/4 ) 
    {
      Location <-factor(myT$Location)
      
      myLm <- lm(myT[,i] ~ Location)
      myAnova <- anova(myLm)
      pValuesLocation[index] <- myAnova$"Pr(>F)"[1]
      
      Name[index] <- names(myT)[i]
      indexes[index] <- index
      absolute_indexes[index]=i
      
      index <- index + 1
    }
  }
  # Making a dataframe out of all the p values and indexes
  dFrame <- data.frame(Name, indexes, absolute_indexes,pValuesLocation)
  # Performing BH adjustments on p values
  dFrame$pValuesLocationAdjusted<- p.adjust( dFrame$pValuesLocation, method = "BH")
  
  dFrame$Name=gsub(pattern ="D_[0-9]__",x = dFrame$Name,replacement = "_")
  dFrame$Name=gsub(pattern ="_Bacteria",x = dFrame$Name,replacement = "Bacteria")
  string <- strsplit(as.character(dFrame$Name),split = "._")
  
  temp_string=do.call(rbind,string)
  dFrame <- cbind(temp_string,dFrame)
  
  colnames(dFrame)[colnames(dFrame)=="1"] <- "Kingdom"
  colnames(dFrame)[colnames(dFrame)=="2"] <- "Phylum"
  colnames(dFrame)[colnames(dFrame)=="3"] <- "Class"
  colnames(dFrame)[colnames(dFrame)=="4"] <- "Order"
  colnames(dFrame)[colnames(dFrame)=="5"] <- "Family"
  colnames(dFrame)[colnames(dFrame)=="6"] <- "Genus"
  colnames(dFrame)[colnames(dFrame)=="7"] <- "Species"
  
  dir.create(paste0(output, input[z]), showWarnings = FALSE)
  outputDir <- paste0(output, input[z],"/")
  dir.create(paste0(outputDir,"PDFs/"), showWarnings = FALSE)
  PDFdir <- paste0(outputDir, "PDFs/")
  
  # Saving results as a text file
  write.table(dFrame, file=paste0(outputDir, input[z],"_GLM_NormalizedCFU_Culture_",loc1,"v",loc2,"_ALL.tsv"), sep="\t",row.names=FALSE)
  if (min(dFrame$pValuesLocationAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesLocationAdjusted < 0.05,]
    Location <- paste(loc1, loc2, sep = " - ")
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesLocationAdjusted
    dFrame2 <- cbind(Taxa, Location, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ALL.tsv"), sep="\t",row.names=FALSE)
  }
  
}
