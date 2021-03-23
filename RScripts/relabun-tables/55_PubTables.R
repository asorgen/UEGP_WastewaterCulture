#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-06-21
#Description: Linear models on normalized counts comparing culture and environmental samples

library(vegan)

rm(list=ls())

pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
output = file.path(moduleDir,"output/")
taxa <- "Filtered_Family"


##### Antibiotic (Table S10) #####
AntibioticModule = dir(pipeRoot, pattern="Antibiotic_RelAbun", full.names=TRUE)
AntibioticPath = file.path(AntibioticModule,"output/")
myT=read.table(paste0(AntibioticPath, taxa, '_Abundance_by_Antibiotic.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

myT$Antibiotic <- factor(myT$Antibiotic, levels = c("Neg", "Amp", "Cip", "Dox", "Sulf"))
var1 <- unique(myT$Antibiotic)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Antibiotic <- paste0(var1[i])
  df1 <- myT[myT$Antibiotic %in% Antibiotic, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Antibiotic, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Antibiotic, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_AntibioticTable.tsv"), sep="\t",row.names=FALSE)




##### Site (Table S3) #####
SiteModule = dir(pipeRoot, pattern="Site_RelAbun", full.names=TRUE)
SitePath = file.path(SiteModule,"output/")
myT=read.table(paste0(SitePath, taxa, '_Abundance_by_Site.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- unique(myT$Site)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Site <- paste0(var1[i])
  df1 <- myT[myT$Site %in% Site, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Site, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Site, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_SiteTable.tsv"), sep="\t",row.names=FALSE)








##### Media #####
MediaModule = dir(pipeRoot, pattern="Media_RelAbun", full.names=TRUE)
MediaPath = file.path(MediaModule,"output/")
myT=read.table(paste0(MediaPath, taxa, '_Abundance_by_Media.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- unique(myT$Media)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Media <- paste0(var1[i])
  df1 <- myT[myT$Media %in% Media, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Media, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Media, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_MediaTable.tsv"), sep="\t",row.names=FALSE)





##### Temperature #####
TemperatureModule = dir(pipeRoot, pattern="Temperature_RelAbun", full.names=TRUE)
TemperaturePath = file.path(TemperatureModule,"output/")
myT=read.table(paste0(TemperaturePath, taxa, '_Abundance_by_Temperature.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- unique(myT$Temperature)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Temperature <- paste0(var1[i])
  df1 <- myT[myT$Temperature %in% Temperature, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Temperature, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Temperature, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_TemperatureTable.tsv"), sep="\t",row.names=FALSE)





##### Location (Table S6) #####
LocationModule = dir(pipeRoot, pattern="Location_RelAbun", full.names=TRUE)
LocationPath = file.path(LocationModule,"output/")
myT=read.table(paste0(LocationPath, taxa, '_Abundance_by_Location.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family
myT$Location <- factor(myT$Location, levels = c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))

var1 <- unique(myT$Location)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Location <- paste0(var1[i])
  df1 <- myT[myT$Location %in% Location, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Location, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Location, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_LocationTable.tsv"), sep="\t",row.names=FALSE)






##### Location ARB (Table S7) #####
LocationModule = dir(pipeRoot, pattern="Location_RelAbun", full.names=TRUE)
LocationPath = file.path(LocationModule,"output/")
myT=read.table(paste0(LocationPath, taxa, '_ARB_Abundance_by_Location.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family
myT$Location <- factor(myT$Location, levels = c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))

var1 <- unique(myT$Location)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Location <- paste0(var1[i])
  df1 <- myT[myT$Location %in% Location, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Location, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Location, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_ARB_LocationTable.tsv"), sep="\t",row.names=FALSE)






##### Site ARB #####
SiteModule = dir(pipeRoot, pattern="Site_RelAbun", full.names=TRUE)
SitePath = file.path(SiteModule,"output/")
myT=read.table(paste0(SitePath, taxa, '_ARB_Abundance_by_Site.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- unique(myT$Site)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Site <- paste0(var1[i])
  df1 <- myT[myT$Site %in% Site, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Site, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Site, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_ARB_SiteTable.tsv"), sep="\t",row.names=FALSE)






##### Media ARB #####
MediaModule = dir(pipeRoot, pattern="Media_RelAbun", full.names=TRUE)
MediaPath = file.path(MediaModule,"output/")
myT=read.table(paste0(MediaPath, taxa, '_ARB_Abundance_by_Media.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- unique(myT$Media)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Media <- paste0(var1[i])
  df1 <- myT[myT$Media %in% Media, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Media, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Media, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_ARB_MediaTable.tsv"), sep="\t",row.names=FALSE)






##### Temperature ARB #####
TemperatureModule = dir(pipeRoot, pattern="Temperature_RelAbun", full.names=TRUE)
TemperaturePath = file.path(TemperatureModule,"output/")
myT=read.table(paste0(TemperaturePath, taxa, '_ARB_Abundance_by_Temperature.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- unique(myT$Temperature)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Temperature <- paste0(var1[i])
  df1 <- myT[myT$Temperature %in% Temperature, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Temperature, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Temperature, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_ARB_TemperatureTable.tsv"), sep="\t",row.names=FALSE)






##### Condition (Table S16) #####
ConditionModule = dir(pipeRoot, pattern="Condition_RelAbun", full.names=TRUE)
ConditionPath = file.path(ConditionModule,"output/")
myT=read.table(paste0(ConditionPath, taxa, '_Abundance_by_Condition.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- unique(myT$Condition)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Condition <- paste0(var1[i])
  df1 <- myT[myT$Condition %in% Condition, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Condition, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Condition, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_ConditionTable.tsv"), sep="\t",row.names=FALSE)


##### Location_Antibiotic #####
Location_AntibioticnModule = dir(pipeRoot, pattern="RelAbun_Location_Antibiotic", full.names=TRUE)
Location_AntibioticPath = file.path(Location_AntibioticnModule,"output/")
myT=read.table(paste0(Location_AntibioticPath, taxa, '_Abundance_by_Location_Antibiotic.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA")

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Location <- paste0(var1[i])
  df1 <- myT[myT$Location %in% Location, ]
  df1$Antibiotic <- factor(df1$Antibiotic, levels = c("Amp","Cip","Dox","Sulf", "Neg"))
  var2 <- unique(df1$Antibiotic)
  
  
  for (j in 1:length(var2)) {
    Antibiotic <- paste0(var2[j])
    df2 <- df1[df1$Antibiotic %in% Antibiotic,]
    RA <- df2$Percentage
    Avg <- df2$Average
    dFrame <- cbind(dFrame, Avg, RA)
    colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Location, Antibiotic, "RelAbun", sep = "_")
    colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Location, Antibiotic, "Average", sep = "_")
  }
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)


write.table(dFrame, file=paste0(output, taxa,"_Location_AntibioticTable.tsv"), sep="\t",row.names=FALSE)



##### UPA & DSA (Table S8) #####
Location_AntibioticnModule = dir(pipeRoot, pattern="RelAbun_Location_Antibiotic", full.names=TRUE)
Location_AntibioticPath = file.path(Location_AntibioticnModule,"output/")
myT=read.table(paste0(Location_AntibioticPath, taxa, '_Abundance_by_Location_Antibiotic.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- c("UPA","DSA")

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Location <- paste0(var1[i])
  df1 <- myT[myT$Location %in% Location, ]
  df1$Antibiotic <- factor(df1$Antibiotic, levels = c("Amp","Cip","Dox","Sulf", "Neg"))
  var2 <- unique(df1$Antibiotic)
  
  
  for (j in 1:length(var2)) {
    Antibiotic <- paste0(var2[j])
    df2 <- df1[df1$Antibiotic %in% Antibiotic,]
    RA <- df2$Percentage
    Avg <- df2$Average
    dFrame <- cbind(dFrame, Avg, RA)
    colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Location, Antibiotic, "RelAbun", sep = "_")
    colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Location, Antibiotic, "Average", sep = "_")
  }
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)


write.table(dFrame, file=paste0(output, taxa,"_UPAvDSATable.tsv"), sep="\t",row.names=FALSE)



##### RES & HOS (Table S9) #####
Location_AntibioticnModule = dir(pipeRoot, pattern="RelAbun_Location_Antibiotic", full.names=TRUE)
Location_AntibioticPath = file.path(Location_AntibioticnModule,"output/")
myT=read.table(paste0(Location_AntibioticPath, taxa, '_Abundance_by_Location_Antibiotic.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

var1 <- c("RES","HOS")

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  Location <- paste0(var1[i])
  df1 <- myT[myT$Location %in% Location, ]
  df1$Antibiotic <- factor(df1$Antibiotic, levels = c("Amp","Cip","Dox","Sulf", "Neg"))
  var2 <- unique(df1$Antibiotic)
  
  
  for (j in 1:length(var2)) {
    Antibiotic <- paste0(var2[j])
    df2 <- df1[df1$Antibiotic %in% Antibiotic,]
    RA <- df2$Percentage
    Avg <- df2$Average
    dFrame <- cbind(dFrame, Avg, RA)
    colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(Location, Antibiotic, "RelAbun", sep = "_")
    colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(Location, Antibiotic, "Average", sep = "_")
  }
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)


write.table(dFrame, file=paste0(output, taxa,"_RESvHOSTable.tsv"), sep="\t",row.names=FALSE)





##### SampleType (Table S10) #####
SampleTypeModule = dir(pipeRoot, pattern="SampleType_RelAbun", full.names=TRUE)
SampleTypePath = file.path(SampleTypeModule,"output/")
myT=read.table(paste0(SampleTypePath, taxa, '_Abundance_by_SampleType.tsv'),sep='\t',header = TRUE)
myT2 <- myT[!duplicated(myT$Family), ]
Phylum <- myT2$Phylum
Class <- myT2$Class
Family <- myT2$Family

myT$SampleType <- factor(myT$SampleType)
var1 <- unique(myT$SampleType)

dFrame <- data.frame(Phylum, Class, Family)
for (i in 1:length(var1)) {
  SampleType <- paste0(var1[i])
  df1 <- myT[myT$SampleType %in% SampleType, ]
  RA <- df1$Percentage
  Avg <- df1$Average
  dFrame <- cbind(dFrame, Avg, RA)
  colnames(dFrame)[colnames(dFrame)=="RA"] <- paste(SampleType, "RelAbun", sep = "_")
  colnames(dFrame)[colnames(dFrame)=="Avg"] <- paste(SampleType, "Average", sep = "_")
}

df <- dFrame[ , grepl( "RelAbun" , names( dFrame ) ) ]
index <- 1
del.rows <- vector()

for (z in 1:nrow(dFrame)) {
  max.relab <- max(df[z,])
  
  if (max.relab < 0.01) {
    del.rows[index] <- Family[z]
    index <- index + 1
  }
}


df2 <- dFrame[dFrame$Family %in% c(del.rows,"Other") ,]
df3 <- df2[ , grepl( "_" , names( dFrame ) ) ]

relabs <- colSums(df3[,1:ncol(df3)])
Other <- c("Other", "Other", "Other", relabs)

dFrame3 <- dFrame[!(dFrame$Family %in% c(del.rows,"Other") ),]
dFrame <- rbind(dFrame3, Other)



write.table(dFrame, file=paste0(output, taxa,"_SampleTypeTable.tsv"), sep="\t",row.names=FALSE)




