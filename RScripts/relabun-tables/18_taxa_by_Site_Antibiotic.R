#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-07-21
#Description: Generate relative abundance tables at each taxonomic level based on site and antibiotic.

rm(list=ls())
taxa=c("Phylum","Class","Order","Filtered_Family","Genus")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Gather", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

L <- 1


for (z in 1:length(taxa)) {
  
  frame=read.table(paste0(inputPath, taxa[z],'_table.tsv'),sep='\t',header = TRUE)
  frame <- frame[frame$SampleType == "Culture",]
  L <- L + 1
  frame$Taxa <- paste(frame$Kingdom,frame$Phylum,frame$Class,frame$Order,frame$Family,frame$Genus,sep = ">")
  list <- c("Mallard", "Sugar")
  
  dFrame=data.frame()
  for (x in 1:length(list)) {
    Site=paste(list[x])
    df=frame[frame$Site %in% Site,]
    
    list2 = factor(unique(df$Antibiotic))
    
    for (y in 1:length(list2)) {
     
       Antibiotic=paste(list2[y])
      df2=df[df$Antibiotic %in% Antibiotic,]
      
      TAXA=factor(unique(frame$Taxa))
      Total_CFU=sum(df2$Abs_Abun)
      
      for (i in 1:length(TAXA)) {
       
         Taxa=paste(TAXA[i])
        df3=df2[df2$Taxa %in% Taxa,]
        Taxa_Total=sum(df3$Abs_Abun)
        Average=mean(df3$Abs_Abun)
        SD=sd(df3$Abs_Abun)
        Percentage=Taxa_Total / Total_CFU
        addition=data.frame(Site,Taxa,Taxa_Total,Total_CFU,Percentage,Average,SD)
        dFrame=rbind(dFrame,addition)

      }
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
  
  write.table(dFrame, file=paste0(output, taxa[z],"_Abundance_by_Site_Antibiotic.tsv"),sep = "\t", row.names = FALSE)
  
}



taxa=c("Phylum","Class","Order","Family","Genus")
L <- 1

for (z in 1:length(taxa)) {
  
  frame=read.table(paste0(inputPath, taxa[z],'_table.tsv'),sep='\t',header = TRUE)
  frame <- frame[frame$SampleType == "Culture",]
  L <- L + 1
  frame$Taxa <- paste(frame$Kingdom,frame$Phylum,frame$Class,frame$Order,frame$Family,frame$Genus,sep = ">")
  list <- c("Mallard", "Sugar")
  
  dFrame=data.frame()
  for (x in 1:length(list)) {
    Site=paste(list[x])
    df=frame[frame$Site %in% Site,]
    
    list2 = factor(unique(df$Antibiotic))
    
    for (y in 1:length(list2)) {
      
      Antibiotic=paste(list2[y])
      df2=df[df$Antibiotic %in% Antibiotic,]
      
      TAXA=factor(unique(frame$Taxa))
      Total_CFU=sum(df2$Abs_Abun)
      
      for (i in 1:length(TAXA)) {
        
        Taxa=paste(TAXA[i])
        df3=df2[df2$Taxa %in% Taxa,]
        Taxa_Total=sum(df3$Abs_Abun)
        Average=mean(df3$Abs_Abun)
        SD=sd(df3$Abs_Abun)
        Percentage=Taxa_Total / Total_CFU
        addition=data.frame(Site,Taxa,Taxa_Total,Total_CFU,Percentage,Average,SD)
        dFrame=rbind(dFrame,addition)
        
      }
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
  
  write.table(dFrame, file=paste0(output, taxa[z],"_Abundance_by_Site_Antibiotic.tsv"),sep = "\t", row.names = FALSE)
  
}
