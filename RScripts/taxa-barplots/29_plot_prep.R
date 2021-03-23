#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-08-21
#Description: 

rm(list=ls())


## .) Set libraries
library(ggplot2)
# library(grid)
# library(scales)



## .) Set directories
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="Gather", full.names=TRUE)
inputPath = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")




## .) Modify abundance by location and antibiotic
taxa_input=c("Filtered_Family")

for (z in 1:length(taxa_input)) {

  
  frame=read.table(paste0(inputPath, taxa_input[z],'_table.tsv'),sep='\t',header = TRUE)
  frame <- frame[frame$SampleType == "Culture",]
  frame$Taxa <- paste(frame$Kingdom,frame$Phylum,frame$Class,frame$Order,frame$Family,frame$Genus,sep = "//")
  
  list <- c("UPA", "RES", 
            "HOS", "INF", 
            "PCI", "PCE", 
            "ATE", "FCE",
            "UV",  "DSA")
  
  dFrame=data.frame()
  for (x in 1:length(list)) {
    Location=paste(list[x])
    df=frame[frame$Location %in% Location,]
    
    list2 = factor(unique(df$Antibiotic))
    
    for (y in 1:length(list2)) {
      
      Antibiotic=paste(list2[y])
      df2=df[df$Antibiotic %in% Antibiotic,]
      
      TAXA=factor(unique(df2$Taxa))
      Total_CFU=sum(df2$Abs_Abun)
      
      for (i in 1:length(TAXA)) {
        
        Taxa=paste(TAXA[i])
        df3=df2[df2$Taxa %in% Taxa,]
        Taxa_Total=sum(df3$Abs_Abun)
        Average=mean(df3$Abs_Abun)
        SD=sd(df3$Abs_Abun)
        Percentage=Taxa_Total / Total_CFU
        addition=data.frame(Location,Antibiotic, Taxa,Taxa_Total,Total_CFU,Percentage,Average,SD)
        dFrame=rbind(dFrame,addition)
        
      }
    }
  }
  
  # write.table(dFrame, file=paste0(output, taxa_input[z],"_Abundance_by_Location_Antibiotic.tsv"),sep = "\t", row.names = FALSE)
}

taxa_input=c("Filtered_Family")

string <- strsplit(as.character(dFrame$Taxa),split = "//")
temp_string=do.call(rbind,string)

dFrame <- cbind(temp_string,dFrame)

colnames(dFrame)[colnames(dFrame)=="1"] <- "Kingdom"
colnames(dFrame)[colnames(dFrame)=="2"] <- "Phylum"
colnames(dFrame)[colnames(dFrame)=="3"] <- "Class"
colnames(dFrame)[colnames(dFrame)=="4"] <- "Order"
colnames(dFrame)[colnames(dFrame)=="5"] <- "Family"

dFrame$Class_Family = paste(dFrame$Class,dFrame$Family,sep = " - ")
dFrame$Class_Family=gsub("Other - Other","Other",dFrame$Class_Family)
dFrame$Class_Family=gsub("Bacilli - Family.XII","Bacilli - Family XII",dFrame$Class_Family)

write.table(dFrame, file=paste0(output, taxa_input,"_Abundance_by_Location_Antibiotic_Figs.tsv"),
            sep = "\t", row.names = FALSE)





## .) Generate family color scheme

taxaNames <- c("Actinobacteria - Corynebacteriaceae","Actinobacteria - Microbacteriaceae","Actinobacteria - Micrococcaceae","Actinobacteria - Mycobacteriaceae","Actinobacteria - Nocardiaceae","Actinobacteria - Propionibacteriaceae","Actinobacteria - Sanguibacteraceae","Actinobacteria - Sporichthyaceae","Alphaproteobacteria - Acetobacteraceae","Alphaproteobacteria - Bradyrhizobiaceae","Alphaproteobacteria - Brucellaceae","Alphaproteobacteria - Caulobacteraceae","Alphaproteobacteria - Rhodobacteraceae","Alphaproteobacteria - Rhodospirillaceae","Alphaproteobacteria - Sphingomonadaceae","Bacilli - Aerococcaceae","Bacilli - Bacillaceae","Bacilli - Enterococcaceae","Bacilli - Family XII","Bacilli - Paenibacillaceae","Bacilli - Planococcaceae","Bacilli - Streptococcaceae","Bacteroidia - Bacteroidaceae","Bacteroidia - Porphyromonadaceae","Betaproteobacteria - Alcaligenaceae","Betaproteobacteria - Burkholderiaceae","Betaproteobacteria - Comamonadaceae","Betaproteobacteria - Methylophilaceae","Betaproteobacteria - Neisseriaceae","Betaproteobacteria - Oxalobacteraceae","Betaproteobacteria - Rhodocyclaceae","Clostridia - Lachnospiraceae","Clostridia - Peptostreptococcaceae","Clostridia - Ruminococcaceae","Cyanobacteria - Family I","Cytophagia - Cytophagaceae","Deltaproteobacteria - Desulfovibrionaceae","Epsilonproteobacteria - Campylobacteraceae","Epsilonproteobacteria - Helicobacteraceae","Flavobacteriia - Cryomorphaceae","Flavobacteriia - Flavobacteriaceae","Flavobacteriia - NS9 marine group","Fusobacteriia - Leptotrichiaceae","Gammaproteobacteria - Aeromonadaceae","Gammaproteobacteria - Chromatiaceae","Gammaproteobacteria - Enterobacteriaceae","Gammaproteobacteria - Halomonadaceae","Gammaproteobacteria - Legionellaceae","Gammaproteobacteria - Moraxellaceae","Gammaproteobacteria - Oceanospirillaceae","Gammaproteobacteria - Pasteurellaceae","Gammaproteobacteria - Pseudomonadaceae","Gammaproteobacteria - Shewanellaceae","Gammaproteobacteria - Thiotrichaceae","Gammaproteobacteria - Vibrionaceae","Gammaproteobacteria - Xanthomonadaceae","Negativicutes - Veillonellaceae","Nitrospira - Nitrospiraceae","Other","Sphingobacteriia - Chitinophagaceae","Sphingobacteriia - NS11-12 marine group","Sphingobacteriia - Sphingobacteriaceae")

hexID <- c("#D15B4D","#008F9C","#7BC853","#FF74A3","#00FFC6","#BDD394","#B0FF61","#014753","#E95EBE","#FFB167","#968AE9","#ECC3E7","#669E46","#91D0CC","#FFE502","#FFFEC8","#009BFF","#C5A9ED","#9CE48B","#43002C","#02FF78","#01D1FF","#F4B187","#D0EF94","#ED672B","#A7573F","#E06B7B","#E7E6E6","#6697E0","#DFFF73","#FF6E41","#0529FF","#DF957D","#A6CACA","#FFF0DA","#B501FF","#3A3838","#5FAD4E","#BF8F00","#E66FFF","#620D00","#223768","#1C5AFF","#005F38","#9F008F","#B0EF80","#7544B1","#5A4817","#001544","#E8CED1","#BB8800","#903CE0","#0E4CA1","#C00000","#C28C9F","#53B63A","#99FF52","#FF900B","#FBE4C5","#6A6882","#BFD62B","#788230")

colorIDs <- data.frame(taxaNames, hexID)
colorIDs <- colorIDs[order(taxaNames),]
write.table(colorIDs, file=paste0(output, "colorIDs.tsv"),sep = "\t", row.names = FALSE)

