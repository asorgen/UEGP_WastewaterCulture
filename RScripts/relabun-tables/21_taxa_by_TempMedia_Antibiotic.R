#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: 

rm(list=ls())
taxa=c("Phylum","Class","Order","Family","Genus", "Filtered_Family")

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Gather", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

L <- 1

for (z in 1:length(taxa)) {
  frame=read.table(paste(inputPath, taxa[z],'_table.tsv',sep=""),sep='\t',header = TRUE)
  frame <- frame[frame$SampleType == "Culture",]
  Condition <- paste(frame$Temperature,frame$Media,sep = "_")
  frame <- cbind(Condition,frame)
  frame$Condition <- factor(frame$Condition, levels = c("BT_LB","RT_LB","BT_R2A","RT_R2A"))
  frame$Antibiotic <- factor(frame$Antibiotic, levels = c("Amp","Cip","Dox","Sulf","Neg"))
  condition_labs=c("LB agar at \n Body temp","LB agar at \n Room temp","R2A agar at \n Body temp","R2A agar at \n Room temp")
  names(condition_labs) <- c("BT_LB","RT_LB","BT_R2A","RT_R2A")
  antibiotic_labs=c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","No Antibiotic")
  names(antibiotic_labs) <- c("Amp","Cip","Dox","Sulf","Neg")
  Locations <- c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA")
  L <- L + 1
  frame$Taxa <- paste(frame$Kingdom,frame$Phylum,frame$Class,frame$Order,frame$Family,frame$Genus,sep = ">")
  
  if (taxa[z] == "Filtered_Family") {
    L <- 5
  }
  
  
  dFrame=data.frame()
  var1 <- factor(unique(frame$Condition))
    
  for (x in 1:length(var1)) {
    Condition=paste(var1[x])
    df=frame[frame$Condition %in% Condition,]
    var2=factor(unique(df$Antibiotic))
    
    for (y in 1:length(var2)) {
      Antibiotic=paste(var2[y])
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
        addition=data.frame(Condition,Antibiotic,Taxa,Taxa_Total,Total_CFU,Percentage,Average,SD)
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
  
  write.table(dFrame, file=paste0(output, taxa[z],"_Abundance_by_TempMedia_and_Antibiotic.tsv"),sep = "\t", row.names = FALSE)
}
