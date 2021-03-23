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
  myT=read.table(paste0(inputDir, input[z],'_Normalized_CFU.tsv'),sep='\t',header = TRUE)
  myT <- myT[myT$SampleType == "Culture",]
  myT <- myT[myT$Antibiotic %in% c("Amp","Cip","Dox","Sulf"),]
  
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
      Antibiotic <- factor(myT$Antibiotic, levels = c("Amp","Cip","Dox","Sulf"))
      
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
  
   # Saving results as a text file
  write.table(dFrame, file=paste0(outputDir, input[z],"_GLM_NormalizedCFU_Culture_ARB_only.tsv"), sep="\t",row.names=FALSE)
  
  if (min(dFrame$pValuesSiteAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesSiteAdjusted < 0.05,]
    Site <- df2$Site_direction
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesSiteAdjusted
    dFrame2 <- cbind(Taxa, Site, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ARB_Site.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesLocationAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesLocationAdjusted < 0.05,]
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesLocationAdjusted
    dFrame2 <- cbind(Taxa, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ARB_Location.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesTemperatureAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesTemperatureAdjusted < 0.05,]
    Temperature <- df2$Temperature_direction
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesTemperatureAdjusted
    dFrame2 <- cbind(Taxa, Temperature, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ARB_Temperature.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesMediaAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesMediaAdjusted < 0.05,]
    Media <- df2$Media_direction
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesMediaAdjusted
    dFrame2 <- cbind(Taxa, Media, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ARB_Media.tsv"), sep="\t",row.names=FALSE)
  }
  
  if (min(dFrame$pValuesAntibioticAdjusted) < 0.05) {
    df2 <- dFrame[dFrame$pValuesAntibioticAdjusted < 0.05,]
    taxaEnd <- which(colnames(df2) == "Name")-1
    Taxa <- df2[,1:taxaEnd]
    pValue <- df2$pValuesAntibioticAdjusted
    dFrame2 <- cbind(Taxa, pValue)
    write.table(dFrame2, file=paste0(output, input[z],"_Significant_Taxa_ARRB_Antibiotic.tsv"), sep="\t",row.names=FALSE)
  }
  
  
  dir.create(paste0(outputDir,"PDFs/"), showWarnings = FALSE)
  PDFdir <- paste0(outputDir, "PDFs/")

    ########## Boxplot for significance differences in Site
  pdf(paste(PDFdir, input[z],"_Site_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
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
  pdf(paste(PDFdir, input[z],"_Location_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
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
  pdf(paste(PDFdir, input[z],"_Temperature_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
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
  pdf(paste(PDFdir, input[z],"_Media_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
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
  pdf(paste(PDFdir, input[z],"_Antibiotic_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
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
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  

  ########## PCoA for Antibiotic
  color <- c("deeppink", "darkorange", "blue", "green3")
  col2=adjustcolor(color[factor(myT$Antibiotic)], alpha.f = 1)
  
  pdf(paste(PDFdir, input[z],"_PCoA_Antibiotic_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
      par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
      
      pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                         xlab = pcoa_p[1], ylab = pcoa_p[2])
      points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
      for (n in 1:4) {
        ordiellipse(pcoa12, myT$Antibiotic, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                    col = color[n], show.groups = levels(factor(myT$Antibiotic))[n], label = T, 
                    font = 2, cex = 1)
      }
      legend("topright", c("Ampicillin","Ciprofloxacin", "Doxycycline", "Sulfamethoxazole"),
             col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
  
  
  ########## PCoA for Temperature
  color=c("steelblue","tomato")
  col2=adjustcolor(color[factor(myT$Temperature)], alpha.f = 1)
  
  pdf(paste(PDFdir, input[z],"_PCoA_Temperature_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
  
  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                     xlab = pcoa_p[1], ylab = pcoa_p[2])
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
  for (n in 1:2) {
    ordiellipse(pcoa12, myT$Temperature, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT$Temperature))[n], label = T, 
                font = 2, cex = 1)
  }
  legend("topright", c("Body temp","Room temp"),
         col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
  
  ########## PCoA for Site
  color=c("steelblue","tomato")
  col2=adjustcolor(color[factor(myT$Site)], alpha.f = 1)
  
  pdf(paste(PDFdir, input[z],"_PCoA_Site_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
  
  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                     xlab = pcoa_p[1], ylab = pcoa_p[2])
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
  for (n in 1:2) {
    ordiellipse(pcoa12, myT$Site, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT$Site))[n], label = T, 
                font = 2, cex = 1)
  }
  legend("topright", c("Mallard","Sugar"),
         col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
  
  ########## PCoA for Media
  color=c("steelblue","tomato")
  col2=adjustcolor(color[factor(myT$Media)], alpha.f = 1)
  
  pdf(paste(PDFdir, input[z],"_PCoA_Media_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
  
  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                     xlab = pcoa_p[1], ylab = pcoa_p[2])
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
  for (n in 1:2) {
    ordiellipse(pcoa12, myT$Media, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT$Media))[n], label = T, 
                font = 2, cex = 1)
  }
  legend("topright", c("LB","R2A"),
         col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
  
  ########## PCoA for Location
  color=c("steelblue","tomato","deeppink","darkorange","purple","green3","black","yellow","magenta","turquoise")
  col2=adjustcolor(color[factor(myT$Location)], alpha.f = 1)
  
  pdf(paste(PDFdir, input[z],"_PCoA_Location_NormalizedCFU_Culture_ARB_only.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
  
  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                     xlab = pcoa_p[1], ylab = pcoa_p[2])
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
  for (n in 1:10) {
    ordiellipse(pcoa12, myT$Location, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT$Location))[n], label = T, 
                font = 2, cex = 1)
  }
  legend("topright", c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"),
         col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
}
