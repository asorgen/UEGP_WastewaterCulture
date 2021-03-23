#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-06-21
#Description: Linear models on normalized counts comparing culture and environmental samples

library(vegan)

rm(list=ls())
input=c("Filtered_Family","Phylum","Class","Order","Family","Genus", "Species")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Normalize", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

for (z in 1:length(input)) {
  inputDir <- paste0(inputPath, input[z], "/")
  myT=read.table(paste0(inputDir, input[z],'_relabun.tsv'),sep='\t',header = TRUE)
  myT <- myT[myT$Location %in% c("UPA", "RES", "HOS", "INF", "PCI", "PCE", "ATE", "FCE", "UV", "DSA"),]
  

  Condition <- paste(myT$Temperature, myT$Media, sep = "_")
  myT <- cbind(Condition, myT)
  myT$Condition=gsub(pattern ="n/a_n/a",x = myT$Condition,replacement = "Environ")
  
  pValuesSite <- vector()
  pValuesLocation <- vector()
  pValuesSampleType <- vector()
  pValuesCondition <- vector()
  Name <- vector()
  index <- 1
  indexes <- vector()
  absolute_indexes=vector()
  Site_direction=vector()
  SampleType_direction=vector()
  taxaStart <- which(colnames(myT)=="CFU.mL")+1
  
  
  for(i in taxaStart:ncol(myT)){
    if(sum(myT[,i]!=0) > nrow(myT)*0.1 ) 
    {
      Site <- factor(myT$Site)
      Location <-factor(myT$Location)
      SampleType <- factor(myT$SampleType)
      Condition <- factor(myT$Condition)
      
      myLm <- lm(myT[,i] ~ Site+Location+SampleType+Condition)
      myAnova <- anova(myLm)
      pValuesSite[index] <- myAnova$"Pr(>F)"[1]
      pValuesLocation[index] <- myAnova$"Pr(>F)"[2]
      pValuesSampleType[index] <- myAnova$"Pr(>F)"[3]
      pValuesCondition[index] <- myAnova$"Pr(>F)"[4]
      
      Name[index] <- names(myT)[i]
      indexes[index] <- index
      absolute_indexes[index]=i
      
      if(mean(myT[myT$Site=="Mallard",i])>mean(myT[myT$Site=="Sugar",i])){
        Site_direction[index]="Mallard"
      }
      else if(mean(myT[myT$Site=="Mallard",i])<mean(myT[myT$Site=="Sugar",i])){
        Site_direction[index]="Sugar"
      }
      else{
        Site_direction[index]="same"
      }
      
      if(mean(myT[myT$SampleType=="Culture",i])>mean(myT[myT$SampleType=="Environ",i])){
        SampleType_direction[index]="Culture"
      }
      else if(mean(myT[myT$SampleType=="Culture",i])<mean(myT[myT$SampleType=="Environ",i])){
        SampleType_direction[index]="Environ"
      }
      else{
        SampleType_direction[index]="same"
      }
      
      index <- index + 1
    }
  }
  # Making a dataframe out of all the p values and indexes
  dFrame <- data.frame(Name, indexes, absolute_indexes,pValuesSite,Site_direction,pValuesLocation,pValuesSampleType,SampleType_direction,pValuesCondition)
  # Performing BH adjustments on p values
  dFrame$pValuesSiteAdjusted<- p.adjust( dFrame$pValuesSite, method = "BH")
  dFrame$pValuesLocationAdjusted<- p.adjust( dFrame$pValuesLocation, method = "BH")
  dFrame$pValuesSampleTypeAdjusted<- p.adjust( dFrame$pValuesSampleType, method = "BH")
  dFrame$pValuesConditionAdjusted<- p.adjust( dFrame$pValuesCondition, method = "BH")
  
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
  
  # Saving results as a text file
  write.table(dFrame, file=paste0(outputDir, input[z],"_GLM_relabun_Combined.tsv"), sep="\t",row.names=FALSE)
  
  if (min(dFrame$pValuesSampleTypeAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesSampleTypeAdjusted < 0.05,]
    SampleType <- df2$SampleType_direction
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesSampleTypeAdjusted
    dFrame2 <- cbind(Taxa, SampleType, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_SampleType.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesConditionAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesConditionAdjusted < 0.05,]
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesConditionAdjusted
    dFrame2 <- cbind(Taxa, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_Condition.tsv"), sep="\t",row.names=FALSE)
  }
  

  dir.create(paste0(outputDir,"PDFs/"), showWarnings = FALSE)
  PDFdir <- paste0(outputDir, "PDFs/")
  
  ########## Boxplot for significance differences in Site
  pdf(paste(PDFdir, input[z],"_Site_relabun_Combined.pdf",sep = ""))
  for (j in 1:nrow(dFrame)){
    absolute_index=dFrame$absolute_indexes[j]
    index=dFrame$indexes[j]
    if(dFrame$pValuesSiteAdjusted[j]>0.05){next}
    bug <- myT[,absolute_index]
    Site <- myT$Site
    myFrame <- data.frame(bug, Site)
    boxplot(myT[,absolute_index] ~ Site, main = paste(colnames(myT[absolute_index]),"\n  adjusted significance  ",format(dFrame$pValuesSiteAdjusted[j])), cex.main=0.5)
    stripchart(bug ~ Site,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
  }
  hist(pValuesSite)
  
  dev.off()
  
  ########## Boxplot for significant differences in Location
  pdf(paste(PDFdir, input[z],"_Location_relabun_Combined.pdf",sep = ""))
  for (j in 1:nrow(dFrame)){
    absolute_index=dFrame$absolute_indexes[j]
    index=dFrame$indexes[j]
    if(dFrame$pValuesLocationAdjusted[j]>0.05){next}
    bug <- myT[,absolute_index]
    Location <- myT$Location
    myFrame <- data.frame(bug, Location)
    boxplot(myT[,absolute_index] ~ Location, main = paste(colnames(myT[absolute_index]),"\n  adjusted significance  ",format(dFrame$pValuesLocationAdjusted[j])), cex.main=0.5)
    stripchart(bug ~ Location,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
  }
  hist(pValuesLocation)
  
  dev.off()
  
  ########## Boxplot for significance differences in SampleType
  pdf(paste(PDFdir, input[z],"_SampleType_relabun_Combined.pdf",sep = ""))
  for (j in 1:nrow(dFrame)){
    absolute_index=dFrame$absolute_indexes[j]
    index=dFrame$indexes[j]
    if(dFrame$pValuesSampleTypeAdjusted[j]>0.05){next}
    bug <- myT[,absolute_index]
    SampleType <- myT$SampleType
    myFrame <- data.frame(bug, SampleType)
    boxplot(myT[,absolute_index] ~ SampleType, main = paste(colnames(myT[absolute_index]),"\n  adjusted significance  ",format(dFrame$pValuesSampleTypeAdjusted[j])), cex.main=0.5)
    stripchart(bug ~ SampleType,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
  }
  hist(pValuesSampleType)
  
  dev.off()
  
  ########## Boxplot for significance differences in Condition
  pdf(paste(PDFdir, input[z],"_Condition_relabun_Combined.pdf",sep = ""))
  for (j in 1:nrow(dFrame)){
    absolute_index=dFrame$absolute_indexes[j]
    index=dFrame$indexes[j]
    if(dFrame$pValuesConditionAdjusted[j]>0.05){next}
    bug <- myT[,absolute_index]
    Condition <- myT$Condition
    myFrame <- data.frame(bug, Condition)
    boxplot(myT[,absolute_index] ~ Condition, main = paste(colnames(myT[absolute_index]),"\n  adjusted significance  ",format(dFrame$pValuesConditionAdjusted[j])), cex.main=0.5)
    stripchart(bug ~ Condition,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
  }
  hist(pValuesCondition)
  
  dev.off()
  
  ########## Make PCoA plots
  pcoa=myT
  myMDS <- capscale(pcoa[,taxaStart:(ncol(pcoa))]~1,distance="bray")
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  ########## PCoA for SampleType
  pdf(paste(PDFdir, input[z],"_PCoA_SampleType_relabun_Combined.pdf",sep = ""))
  for (x in 1:5){
    for (y in 2:4){
      if (x==y){break}
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(myMDS$CA$u[,x], myMDS$CA$u[,y],
           xlab=paste("MDS Axis ", x, ": ", format(percentVariance[x],digits=3), sep=""),
           ylab=paste("MDS Axis ", y, ": ", format(percentVariance[y],digits=3), sep=""),
           main=paste("Bray-Curtis Beta Diversity at the",input[z],"Level",sep = " "), cex=2.0,
           col= ifelse(pcoa$SampleType=="Culture", "red","blue"),
           ,bty='L'
      )
      par(xpd=TRUE)
      legend("bottomright",c("Culture","Environmental"),cex=0.7,col=c("red", "blue"), pch=c(16,16))
    }
  }
  dev.off()
  
  ########## PCoA for Site
  pdf(paste(PDFdir, input[z],"_PCoA_Site_relabun_Combined.pdf",sep = ""))
  for (x in 1:5){
    for (y in 2:4){
      if (x==y){break}
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(myMDS$CA$u[,x], myMDS$CA$u[,y],
           xlab=paste("MDS Axis ", x, ": ", format(percentVariance[x],digits=3), sep=""),
           ylab=paste("MDS Axis ", y, ": ", format(percentVariance[y],digits=3), sep=""),
           main=paste("Bray-Curtis Beta Diversity at the",input[z],"Level",sep = " "), cex=2.0,
           col= ifelse(pcoa$Site=="Mallard", "red","blue"),
           ,bty='L'
      )
      par(xpd=TRUE)
      legend("bottomright",c("Mallard","Sugar"),cex=0.7,col=c("red", "blue"), pch=c(16,16))
    }
  }
  dev.off()
  
  
  ########## PCoA for Location
  pdf(paste(PDFdir, input[z],"_PCoA_Location_relabun_Combined.pdf",sep = ""))
  for (x in 1:5){
    for (y in 2:4){
      if (x==y){break}
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(myMDS$CA$u[,x], myMDS$CA$u[,y],
           xlab=paste("MDS Axis ", x, ": ", format(percentVariance[x],digits=3), sep=""),
           ylab=paste("MDS Axis ", y, ": ", format(percentVariance[y],digits=3), sep=""),
           main=paste("Bray-Curtis Beta Diversity at the",input[z],"Level",sep = " "), cex=2.0,
           col= ifelse(pcoa$Location=="UPA", "darkblue",ifelse(pcoa$Location=="RES", "red",ifelse(pcoa$Location=="HOS", "orange",ifelse(pcoa$Location=="INF", "green3",ifelse(pcoa$Location=="PCI", "darkviolet",ifelse(pcoa$Location=="PCE", "mediumpurple",ifelse(pcoa$Location=="ATE", "maroon1",ifelse(pcoa$Location=="FCE", "goldenrod1",ifelse(pcoa$Location=="UV", "darkturquoise","cadetblue1"))))))))),
           ,bty='L'
      )
      par(xpd=TRUE)
      legend("bottomright",c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"),cex=0.7,col=c("darkblue","red","orange","green3","darkviolet","mediumpurple","maroon1","goldenrod1","darkturquoise","aquamarine"), pch=c(16,16,16,16,16,16,16,16,16,16))
    }
  }
  dev.off()
  
  
  ########## PCoA for Condition
  pdf(paste(PDFdir, input[z],"_PCoA_Condition_relabun_Combined.pdf",sep = ""))
  for (x in 1:5){
    for (y in 2:4){
      if (x==y){break}
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(myMDS$CA$u[,x], myMDS$CA$u[,y],
           xlab=paste("MDS Axis ", x, ": ", format(percentVariance[x],digits=3), sep=""),
           ylab=paste("MDS Axis ", y, ": ", format(percentVariance[y],digits=3), sep=""),
           main=paste("Bray-Curtis Beta Diversity at the",input[z],"Level",sep = " "), cex=2.0,
           col= ifelse(pcoa$Condition=="BT_LB", "darkblue",ifelse(pcoa$Condition=="BT_R2A", "red",ifelse(pcoa$Condition=="Environ_", "orange",ifelse(pcoa$Condition=="RT_LB", "green3","darkviolet")))),
           ,bty='L'
      )
      par(xpd=TRUE)
      legend("bottomright",c("BT_LB","BT_R2A","Environ_","RT_LB","RT_R2A"),cex=0.7,col=c("darkblue","red","orange","green3","darkviolet"), pch=c(16,16,16,16,16,16,16,16,16,16))
    }
  }
  dev.off()
  
}




