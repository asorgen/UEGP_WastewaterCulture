#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-08-21
#Description: 

rm(list=ls())


## .) Set libraries
library(vegan)



## .) Set directories
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Normalize", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")




## .) Create PCoAs for each taxonomic level
taxa_input=c("Phylum","Class","Order","Filtered_Family","Family","Genus")

for (g in 1:length(taxa_input)) {
  inputDir <- paste0(inputPath, taxa_input[g],"/",taxa_input[g])
  outputDir <- paste0(output, taxa_input[g])
  loc <- c("UPA","DSA")
  
  myT=read.table(paste(inputDir,"_Normalized_CFU.tsv",sep=""),sep='\t',header = TRUE)
  myT <- myT[myT$SampleType == "Culture",]
  myT <- myT[(myT$Location %in% loc),]
  myT$Location=gsub("DSA","Downstream",myT$Location)
  myT$Location=gsub("UPA","Upstream",myT$Location)
  taxaStart <- which(colnames(myT)=="CFU.mL")+1
  
  
  myT[, taxaStart:ncol(myT)][is.na(myT[, taxaStart:ncol(myT)])] <- 0

  myMDS <- capscale(myT[,taxaStart:(ncol(myT))]~1,distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")
  
  adon.results<-adonis(myT[, taxaStart:ncol(myT)] ~ myT$Location, method="bray",perm=999)
  # print(adon.results)
  capture.output(adon.results, file = paste0(outputDir,"_adonis_UPAvHOS.txt"))
  Title <- paste0(taxa_input[g], " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])
  
  color <- c("steelblue", "tomato")
  col2=adjustcolor(color[factor(myT$Location)], alpha.f = 1)
  
  pdf(paste(outputDir, "_PCoA_UPAvDSA_Normalized_CFU_Culture_ALL.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))
  
  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2, 
                     xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
  for (n in 1:2) {
    ordiellipse(pcoa12, myT$Location, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT$Location))[n], label = T, 
                font = 2, cex = 1)
  }
  legend("topright", c("Upstream","Downstream"),
         col = color, cex = 1.5, pch = 16, bty = "n")
  dev.off()
  
  
}

