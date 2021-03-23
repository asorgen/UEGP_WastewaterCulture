#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-08-21
#Description: Generate PCoAs based on Bray-Curtis distances for OTUs classfied from the genus to phylum levels clustered by sample type.

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

  myT=read.table(paste(inputDir,"_Normalized.tsv",sep=""),sep='\t',header = TRUE)
  taxaStart <- which(colnames(myT)=="CFU.mL")+1

  myT[, taxaStart:ncol(myT)][is.na(myT[, taxaStart:ncol(myT)])] <- 0
  myT$SampleType=gsub("Environ","Environmental",myT$SampleType)

  myMDS <- capscale(myT[,taxaStart:(ncol(myT))]~1,distance="bray")
  percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
  pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")

  color <- c("steelblue", "tomato")
  col2=adjustcolor(color[factor(myT$SampleType)], alpha.f = 1)

  adon.results<-adonis(myT[, taxaStart:ncol(myT)] ~ myT$SampleType, method="bray",perm=999)
  # print(adon.results)
  capture.output(adon.results, file = paste0(outputDir,"_adonis_SampleType.txt"))
  Title <- paste0(taxa_input[g], " PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])

  pdf(paste(outputDir, "_PCoA_SampleType_Normalized.pdf",sep = ""))
  par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))

  pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 1.5,
                     xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
  points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1)
  for (n in 1:2) {
    ordiellipse(pcoa12, myT$SampleType, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
                col = color[n], show.groups = levels(factor(myT$SampleType))[n], label = T,
                font = 2, cex = 1.25)
  }
  legend("topright", c("Culture","Environmental"),
         col = color, cex = 1.25, pch = 16, bty = "n")
  dev.off()


}

inputModule = dir(pipeRoot, pattern="OTUtable", full.names=TRUE)
inputPath = file.path(inputModule,"output/")

myT=read.table(paste(inputPath,"Normalized_OTUtable.tsv",sep=""),sep='\t',header = TRUE)
taxaStart <- which(colnames(myT)=="SequenceFrequency")+1

myT[, taxaStart:ncol(myT)][is.na(myT[, taxaStart:ncol(myT)])] <- 0
myT$SampleType=gsub("Environ","Environmental",myT$SampleType)

myMDS <- capscale(myT[,taxaStart:(ncol(myT))]~1,distance="bray")
percentVariance <- round(eigenvals(myMDS)/sum(eigenvals(myMDS))*100,2)
pcoa_p=paste("PCoA",c(1:5)," (",percentVariance,"%)",sep="")

color <- c("steelblue", "tomato")
col2=adjustcolor(color[factor(myT$SampleType)], alpha.f = 1)

adon.results<-adonis(myT[, taxaStart:ncol(myT)] ~ myT$SampleType, method="bray",perm=999)
# print(adon.results)
capture.output(adon.results, file = paste0(outputDir,"_adonis_SampleType_OTU.txt"))
Title <- paste0("OTU PERMANOVA p = ", adon.results$aov.tab$`Pr(>F)`[1])

pdf(paste(output, "OTU_PCoA_SampleType_Normalized.pdf",sep = ""))
par(mar=c(5, 5, 5, 5), xpd=TRUE, mfrow = c(1,1))

pcoa12 <- ordiplot(myMDS, choices = c(1,2), display = "sites", type = "none", cex.lab = 2,
                   xlab = pcoa_p[1], ylab = pcoa_p[2], main = Title)
points(pcoa12, "sites", col = adjustcolor(col2, alpha.f = 0.5), pch = 16, cex = 1.5)
for (n in 1:2) {
  ordiellipse(pcoa12, myT$SampleType, kind = "se", conf = 0.95, lwd = 4, draw = "lines",
              col = color[n], show.groups = levels(factor(myT$SampleType))[n], label = T,
              font = 2, cex = 1)
}
legend("bottomright", c("Culture","Environmental"),
       col = color, cex = 1.5, pch = 16, bty = "n")
dev.off()

