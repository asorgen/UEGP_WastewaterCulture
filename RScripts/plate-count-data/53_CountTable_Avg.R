#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: 

## Libraries
library(tidyr)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))

moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")

frame=read.table(paste0(input, "counts.txt"),sep="\t",header = TRUE)
frame=frame[frame$TimePoint=="4",]

##### Location Antibiotic Count Table (Table S4) #####
var1 <- unique(frame$Location)
dFrame <- data.frame()
Location <- vector()
Antibiotic <- vector()
Average <- vector()
index <- 1

for (i in 1:length(var1)) {
  Location[index] <- paste0(var1[i])
  df1 <- frame[frame$Location %in% Location, ]
  df1$Antibiotic <- factor(df1$Antibiotic, levels = c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole", "Negative"))
  var2 <- unique(df1$Antibiotic)

  for (j in 1:length(var2)) {
    Antibiotic[index] <- paste0(var2[j])
    df2 <- df1[df1$Antibiotic %in% Antibiotic,]
    Average[index] <- mean(df2$CFU.mL)
    row <- data.frame(Location, Antibiotic, Average)
    dFrame <- rbind(dFrame, row)
    
  }
}

dFrame2 <- spread(dFrame, Antibiotic, Average)



frame2 <- frame[!(frame$Antibiotic %in% "Negative"),]

dFrame <- data.frame()
Location <- vector()
Average <- vector()
index <- 1

for (i in 1:length(var1)) {
  Location[index] <- paste0(var1[i])
  df1 <- frame[frame$Location %in% Location, ]
  Average[index] <- mean(df1$CFU.mL)
  row <- data.frame(Location, Average)
  dFrame <- rbind(dFrame, row)
  
}
colnames(dFrame)[colnames(dFrame)=="Average"] <- "Location Avg"

dFrame2 <- merge(dFrame2, dFrame, by = "Location")

dFrame <- data.frame()
Location <- vector()
Average <- vector()
index <- 1

for (i in 1:length(var1)) {
  Location[index] <- paste0(var1[i])
  df1 <- frame2[frame2$Location %in% Location, ]
  Average[index] <- mean(df1$CFU.mL)
  row <- data.frame(Location, Average)
  dFrame <- rbind(dFrame, row)
  
}
colnames(dFrame)[colnames(dFrame)=="Average"] <- "Location Avg (ARB only)"

dFrame2 <- merge(dFrame2, dFrame, by = "Location")
Antibiotics <- c("Antibiotic Avg", mean(dFrame2$Ampicillin), mean(dFrame2$Ciprofloxacin), mean(dFrame2$Doxycycline), mean(dFrame2$Negative), mean(dFrame2$Sulfamethoxazole), mean(dFrame2$ALL), mean(dFrame2$ARB))

final <- rbind(dFrame2, Antibiotics)
write.table(final, file=paste0(output, "Location_Antibiotic_CountTable.tsv"), sep="\t",row.names=FALSE)







##### Temperature Media Count Table (Table S11) #####
var1 <- unique(frame$Temperature)
dFrame <- data.frame()
Temperature <- vector()
Media <- vector()
Average <- vector()
index <- 1

for (i in 1:length(var1)) {
  Temperature[index] <- paste0(var1[i])
  df1 <- frame[frame$Temperature %in% Temperature, ]
  df1$Media <- factor(df1$Media)
  var2 <- unique(df1$Media)
  
  for (j in 1:length(var2)) {
    Media[index] <- paste0(var2[j])
    df2 <- df1[df1$Media %in% Media,]
    Average[index] <- mean(df2$CFU.mL)
    row <- data.frame(Temperature, Media, Average)
    dFrame <- rbind(dFrame, row)
    
  }
}

dFrame2 <- spread(dFrame, Media, Average)


dFrame <- data.frame()
Temperature <- vector()
Average <- vector()
index <- 1

for (i in 1:length(var1)) {
  Temperature[index] <- paste0(var1[i])
  df1 <- frame[frame$Temperature %in% Temperature, ]
  Average[index] <- mean(df1$CFU.mL)
  row <- data.frame(Temperature, Average)
  dFrame <- rbind(dFrame, row)
  
}
colnames(dFrame)[colnames(dFrame)=="Average"] <- "Temperature Avg"

dFrame2 <- merge(dFrame2, dFrame, by = "Temperature")


Medias <- c("Media Avg", mean(dFrame2$LB), mean(dFrame2$R2A))

final <- rbind(dFrame2, Medias)
write.table(final, file=paste0(output, "Temperature_Media_CountTable_ALL.tsv"), sep="\t",row.names=FALSE)





var1 <- unique(frame2$Temperature)
dFrame <- data.frame()
Temperature <- vector()
Media <- vector()
Average <- vector()
index <- 1

for (i in 1:length(var1)) {
  Temperature[index] <- paste0(var1[i])
  df1 <- frame2[frame2$Temperature %in% Temperature, ]
  df1$Media <- factor(df1$Media)
  var2 <- unique(df1$Media)
  
  for (j in 1:length(var2)) {
    Media[index] <- paste0(var2[j])
    df2 <- df1[df1$Media %in% Media,]
    Average[index] <- mean(df2$CFU.mL)
    row <- data.frame(Temperature, Media, Average)
    dFrame <- rbind(dFrame, row)
    
  }
}

dFrame2 <- spread(dFrame, Media, Average)


dFrame <- data.frame()
Temperature <- vector()
Average <- vector()
index <- 1

for (i in 1:length(var1)) {
  Temperature[index] <- paste0(var1[i])
  df1 <- frame2[frame2$Temperature %in% Temperature, ]
  Average[index] <- mean(df1$CFU.mL)
  row <- data.frame(Temperature, Average)
  dFrame <- rbind(dFrame, row)
  
}
colnames(dFrame)[colnames(dFrame)=="Average"] <- "Temperature Avg"

dFrame2 <- merge(dFrame2, dFrame, by = "Temperature")


Medias <- c("Media Avg", mean(dFrame2$LB), mean(dFrame2$R2A))

final <- rbind(dFrame2, Medias)
write.table(final, file=paste0(output, "Temperature_Media_CountTable_ARB.tsv"), sep="\t",row.names=FALSE)