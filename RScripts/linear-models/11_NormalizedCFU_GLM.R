#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-06-21
#Description: 

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
  myT=read.table(paste0(inputDir, input[z],'_Normalized_CFU.tsv'),sep='\t',header = TRUE)
  myT <- myT[myT$SampleType == "Culture",]

  pValuesSite <- vector()
  pValuesLocation <- vector()
  pValuesTemperature <- vector()
  pValuesMedia <- vector()
  pValuesAntibiotic <- vector()
  Name <- vector()
  index <- 1
  indexes <- vector()
  absolute_indexes=vector()
  Site_direction=vector()
  Temperature_direction=vector()
  Media_direction=vector()
  taxaStart <- which(colnames(myT)=="CFU.mL")+1
  
  
  for(i in taxaStart:ncol(myT)){
    if(sum(myT[,i]!=0) > nrow(myT)*0.1 ) 
    {
      Site <- factor(myT$Site)
      Location <-factor(myT$Location)
      Temperature <- factor(myT$Temperature)
      Media <- factor(myT$Media)
      Antibiotic <- factor(myT$Antibiotic, levels = c("Neg","Amp","Cip","Dox","Sulf"))
      
      myLm <- lm(myT[,i] ~ Site+Location+Temperature+Media+Antibiotic)
      myAnova <- anova(myLm)
      
      pValuesSite[index] <- myAnova$"Pr(>F)"[1]
      pValuesLocation[index] <- myAnova$"Pr(>F)"[2]
      pValuesTemperature[index] <- myAnova$"Pr(>F)"[3]
      pValuesMedia[index] <- myAnova$"Pr(>F)"[4]
      pValuesAntibiotic[index] <- myAnova$"Pr(>F)"[5]
      
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
      
      if(mean(myT[myT$Temperature=="BT",i])>mean(myT[myT$Temperature=="RT",i])){
        Temperature_direction[index]="BT"
      }
      else if(mean(myT[myT$Temperature=="BT",i])<mean(myT[myT$Temperature=="RT",i])){
        Temperature_direction[index]="RT"
      }
      else{
        Temperature_direction[index]="same"
      }
      
      if(mean(myT[myT$Media=="LB",i])>mean(myT[myT$Media=="R2A",i])){
        Media_direction[index]="LB"
      }
      else if(mean(myT[myT$Media=="LB",i])<mean(myT[myT$Media=="R2A",i])){
        Media_direction[index]="R2A"
      }
      else{
        Media_direction[index]="same"
      }
      
      index <- index + 1
    }
  }
  # Making a dataframe out of all the p values and indexes
  dFrame <- data.frame(Name, indexes, absolute_indexes,pValuesSite,Site_direction,pValuesLocation,pValuesTemperature,Temperature_direction,pValuesMedia,Media_direction,pValuesAntibiotic)
  # Performing BH adjustments on p values
  dFrame$pValuesSiteAdjusted<- p.adjust( dFrame$pValuesSite, method = "BH")
  dFrame$pValuesLocationAdjusted<- p.adjust( dFrame$pValuesLocation, method = "BH")
  dFrame$pValuesTemperatureAdjusted<- p.adjust( dFrame$pValuesTemperature, method = "BH")
  dFrame$pValuesMediaAdjusted<- p.adjust( dFrame$pValuesMedia, method = "BH")
  dFrame$pValuesAntibioticAdjusted<- p.adjust( dFrame$pValuesAntibiotic, method = "BH")
  
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
  write.table(dFrame, file=paste0(outputDir, input[z],"_GLM_NormalizedCFU_Culture.tsv"), sep="\t",row.names=FALSE)
  
  if (min(dFrame$pValuesSiteAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesSiteAdjusted < 0.05,]
    Site <- df2$Site_direction
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesSiteAdjusted
    dFrame2 <- cbind(Taxa, Site, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ALL_Site.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesLocationAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesLocationAdjusted < 0.05,]
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesLocationAdjusted
    dFrame2 <- cbind(Taxa, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ALL_Location.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesTemperatureAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesTemperatureAdjusted < 0.05,]
    Temperature <- df2$Temperature_direction
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesTemperatureAdjusted
    dFrame2 <- cbind(Taxa, Temperature, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ALL_Temperature.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesMediaAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesMediaAdjusted < 0.05,]
    Media <- df2$Media_direction
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesMediaAdjusted
    dFrame2 <- cbind(Taxa, Media, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ALL_Media.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesAntibioticAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesAntibioticAdjusted < 0.05,]
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesAntibioticAdjusted
    dFrame2 <- cbind(Taxa, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ALL_Antibiotic.tsv"), sep="\t",row.names=FALSE)
  }
  
    ########## Boxplot for significance differences in Site
  pdf(paste(PDFdir, input[z],"_Site_NormalizedCFU_Culture.pdf",sep = ""))
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
  
  ########## Boxplot for significance differences in Location
  pdf(paste(PDFdir, input[z],"_Location_NormalizedCFU_Culture.pdf",sep = ""))
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
  
  ########## Boxplot for significance differences in Temperature
  pdf(paste(PDFdir, input[z],"_Temperature_NormalizedCFU_Culture.pdf",sep = ""))
  for (j in 1:nrow(dFrame)){
    absolute_index=dFrame$absolute_indexes[j]
    index=dFrame$indexes[j]
    if(dFrame$pValuesTemperatureAdjusted[j]>0.05){next}
    bug <- myT[,absolute_index]
    Temperature <- myT$Temperature
    myFrame <- data.frame(bug, Temperature)
    boxplot(myT[,absolute_index] ~ Temperature, main = paste(colnames(myT[absolute_index]),"\n adjusted significance  ",format(dFrame$pValuesTemperatureAdjusted[j])), cex.main=0.5)
    stripchart(bug ~ Temperature,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
  }
  hist(pValuesTemperature)
  dev.off()
  
  ########## Boxplot for significance differences in Media
  pdf(paste(PDFdir, input[z],"_Media_NormalizedCFU_Culture.pdf",sep = ""))
  for (j in 1:nrow(dFrame)){
    absolute_index=dFrame$absolute_indexes[j]
    index=dFrame$indexes[j]
    if(dFrame$pValuesMediaAdjusted[j]>0.05){next}
    bug <- myT[,absolute_index]
    Media <- myT$Media
    myFrame <- data.frame(bug, Media)
    boxplot(myT[,absolute_index] ~ Media, main = paste(colnames(myT[absolute_index]),"\n adjusted significance  ",format(dFrame$pValuesMediaAdjusted[j])), cex.main=0.5)
    stripchart(bug ~ Media,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
  }
  hist(pValuesMedia)
  
  dev.off()
  
  ########## Boxplot for significance differences in Antibiotic
  pdf(paste(PDFdir, input[z],"_Antibiotic_NormalizedCFU_Culture.pdf",sep = ""))
  for (j in 1:nrow(dFrame)){
    absolute_index=dFrame$absolute_indexes[j]
    index=dFrame$indexes[j]
    if(dFrame$pValuesAntibioticAdjusted[j]>0.05){next}
    bug <- myT[,absolute_index]
    Antibiotic <- myT$Antibiotic
    myFrame <- data.frame(bug, Antibiotic)
    boxplot(myT[,absolute_index] ~ Antibiotic, main = paste(colnames(myT[absolute_index]),"\n adjusted significance  ",format(dFrame$pValuesAntibioticAdjusted[j])), cex.main=0.5, cex.lab = 0.5)
    stripchart(bug ~ Antibiotic,method="jitter",data = myFrame,vertical = TRUE, pch = 1, add=TRUE )
  }
  hist(pValuesAntibiotic)
  
  dev.off()
  
  
  
  ########## Make PCoA plots
  pcoa=myT
  myMDS <- capscale(pcoa[,taxaStart:(ncol(pcoa))]~1,distance="bray")
  percentVariance <- eigenvals(myMDS)/sum(eigenvals(myMDS))
  
  ########## PCoA for Antibiotic
  pdf(paste(PDFdir, input[z],"_PCoA_Antibiotic_NormalizedCFU_Culture.pdf",sep = ""))
  for (x in 1:5){
    for (y in 2:4){
      if (x==y){break}
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(myMDS$CA$u[,x], myMDS$CA$u[,y],
           xlab=paste("MDS Axis ", x, ": ", format(percentVariance[x],digits=3), sep=""),
           ylab=paste("MDS Axis ", y, ": ", format(percentVariance[y],digits=3), sep=""),
           main=paste("Bray-Curtis Beta Diversity at the",input[z],"Level",sep = " "), cex=2.0,
           col= ifelse(pcoa$Antibiotic=="Amp", "deeppink",ifelse(pcoa$Antibiotic=="Cip", "darkorange",ifelse(pcoa$Antibiotic=="Dox", "blue",ifelse(pcoa$Antibiotic=="Sulf", "green3","black")))),
           ,bty='L'
      )
      par(xpd=TRUE)
      legend("bottomright",c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","Control"),cex=0.7,col=c("deeppink", "darkorange","blue","green3","black"), pch=c(17,16,16,16,16))
    }
  }
  dev.off()
  
  ########## PCoA for Temperature
  pdf(paste(PDFdir, input[z],"_PCoA_Temperature_NormalizedCFU_Culture.pdf",sep = ""))
  for (x in 1:5){
    for (y in 2:4){
      if (x==y){break}
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(myMDS$CA$u[,x], myMDS$CA$u[,y],
           xlab=paste("MDS Axis ", x, ": ", format(percentVariance[x],digits=3), sep=""),
           ylab=paste("MDS Axis ", y, ": ", format(percentVariance[y],digits=3), sep=""),
           main=paste("Bray-Curtis Beta Diversity at the",input[z],"Level",sep = " "), cex=2.0,
           col= ifelse(pcoa$Temperature=="BT", "red","blue"),
           ,bty='L'
      )
      par(xpd=TRUE)
      legend("bottomright",c("Body temp","Room temp"),cex=0.7,col=c("red", "blue"), pch=c(16,16))
    }
  }
  dev.off()
  
  ########## PCoA for Site
  pdf(paste(PDFdir, input[z],"_PCoA_Site_NormalizedCFU_Culture.pdf",sep = ""))
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
  
  ########## PCoA for Media
  pdf(paste(PDFdir, input[z],"_PCoA_Media_NormalizedCFU_Culture.pdf",sep = ""))
  for (x in 1:5){
    for (y in 2:4){
      if (x==y){break}
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(myMDS$CA$u[,x], myMDS$CA$u[,y],
           xlab=paste("MDS Axis ", x, ": ", format(percentVariance[x],digits=3), sep=""),
           ylab=paste("MDS Axis ", y, ": ", format(percentVariance[y],digits=3), sep=""),
           main=paste("Bray-Curtis Beta Diversity at the",input[z],"Level",sep = " "), cex=2.0,
           col= ifelse(pcoa$Media=="LB", "purple","green"),
           ,bty='L'
      )
      par(xpd=TRUE)
      legend("bottomright",c("LB","R2A"),cex=0.7,col=c("purple", "green"), pch=c(16,16))
    }
  }
  dev.off()
  
  ########## PCoA for Location
  pdf(paste(PDFdir, input[z],"_PCoA_Location_NormalizedCFU_Culture.pdf",sep = ""))
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
  
}


