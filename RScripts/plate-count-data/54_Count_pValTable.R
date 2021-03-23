#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: Generate table containing all p values for colony count analyses.

## Libraries
library(tidyr)
library(data.table)
library(dplyr)

rm(list=ls())
##### Prep #####
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")

frame=read.table(paste0(input, "counts.txt"),sep="\t",header = TRUE)
frame=frame[frame$TimePoint=="4",]
frame$CFU <- frame$CFU.mL + 1
frame$logCFU_mL <- log10(frame$CFU)

frame$Site <- factor(frame$Site)
frame$Location <- factor(frame$Location)
frame$Temperature <- factor(frame$Temperature)
frame$Media <- factor(frame$Media)
frame$Antibiotic <- factor(frame$Antibiotic)

##### ARB p values #####
frame2 <- frame[!(frame$Antibiotic == "Negative"),]
frame2$Site <- factor(frame2$Site)
frame2$Location <- factor(frame2$Location)
frame2$Temperature <- factor(frame2$Temperature)
frame2$Media <- factor(frame2$Media)
frame2$Antibiotic <- factor(frame2$Antibiotic)


my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Site","Mallard - ",df$terms)
df$terms = gsub("Location","ATE - ",df$terms)


## Site
df2 <- df[ df$terms %like% "Mallard" ,]
Terms <- df2$terms
ARB <- df2$Pr...t..
final2 <- data.frame(Terms, ARB)


## Location
var2 <- c("Influent", "Primary Clarifier Influent", "Primary Clarifier Effluent", "Final Clarifier Effluent", "UV Treated")
for (i in 1:length(var2)) {
  LocTerms <- paste0("ATE - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "UV Treated")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","UV - ",df$terms)
var2 <- c( "Final Clarifier Effluent", "Downstream")
for (i in 1:length(var2)) {
  LocTerms <- paste0("UV - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "Influent")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","INF - ",df$terms)
var2 <- c( "Residential Sewage", "Hospital Sewage")
for (i in 1:length(var2)) {
  LocTerms <- paste0("INF - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "Downstream")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","DSA - ",df$terms)
var2 <- c( "Upstream")
for (i in 1:length(var2)) {
  LocTerms <- paste0("DSA - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "Residential Sewage")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","RES - ",df$terms)
var2 <- c( "Hospital Sewage")
for (i in 1:length(var2)) {
  LocTerms <- paste0("RES - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "Primary Clarifier Influent")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","PCI - ",df$terms)
var2 <- c( "Primary Clarifier Effluent")
for (i in 1:length(var2)) {
  LocTerms <- paste0("PCI - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}


## Temperature
df$terms = gsub("Temperature","BT - ",df$terms)
df2 <- df[ df$terms %like% "BT" ,]
Terms <- df2$terms
ARB <- df2$Pr...t..
final2 <- rbind(final2, c(Terms, ARB))


## Media
df$terms = gsub("Media","LB - ",df$terms)
df2 <- df[ df$terms %like% "LB" ,]
Terms <- df2$terms
ARB <- df2$Pr...t..
final <- rbind(final2, c(Terms, ARB))


##### ARB p values #####
var1 <- c("Negative", "Ampicillin", "Ciprofloxacin", "Doxycycline", "Sulfamethoxazole")

for (j in 1:length(var1)) {
  Antibiotic <- paste0(var1[j])
  frame2 <- frame[frame$Antibiotic %in% Antibiotic,]
  frame2$Site <- factor(frame2$Site)
  frame2$Location <- factor(frame2$Location)
  frame2$Temperature <- factor(frame2$Temperature)
  frame2$Media <- factor(frame2$Media)

  
  my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media, data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Site","Mallard - ",df$terms)
  df$terms = gsub("Location","ATE - ",df$terms)
  
  
  ## Site
  df2 <- df[ df$terms %like% "Mallard" ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- data.frame(Terms, ARB)
  
  
  ## Location
  var2 <- c("Influent", "Primary Clarifier Influent", "Primary Clarifier Effluent", "Final Clarifier Effluent", "UV Treated")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("ATE - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  frame2$Location <- relevel(frame2$Location, "UV Treated")
  my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","UV - ",df$terms)
  var2 <- c( "Final Clarifier Effluent", "Downstream")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("UV - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  frame2$Location <- relevel(frame2$Location, "Influent")
  my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","INF - ",df$terms)
  var2 <- c( "Residential Sewage", "Hospital Sewage")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("INF - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  frame2$Location <- relevel(frame2$Location, "Downstream")
  my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","DSA - ",df$terms)
  var2 <- c( "Upstream")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("DSA - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  frame2$Location <- relevel(frame2$Location, "Residential Sewage")
  my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","RES - ",df$terms)
  var2 <- c( "Hospital Sewage")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("RES - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  frame2$Location <- relevel(frame2$Location, "Primary Clarifier Influent")
  my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","PCI - ",df$terms)
  var2 <- c( "Primary Clarifier Effluent")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("PCI - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  
  ## Temperature
  df$terms = gsub("Temperature","BT - ",df$terms)
  df2 <- df[ df$terms %like% "BT" ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
  
  ## Media
  df$terms = gsub("Media","LB - ",df$terms)
  df2 <- df[ df$terms %like% "LB" ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  colnames(final2)[colnames(final2)=="ARB"] <- Antibiotic[i]
  
  final <- merge(final, final2, by = "Terms")
  
}



##### ALL p values #####
frame$Site <- factor(frame$Site)
frame$Location <- factor(frame$Location)
frame$Temperature <- factor(frame$Temperature)
frame$Media <- factor(frame$Media)
frame$Antibiotic <- factor(frame$Antibiotic)


my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Site","Mallard - ",df$terms)
df$terms = gsub("Location","ATE - ",df$terms)


## Site
df2 <- df[ df$terms %like% "Mallard" ,]
Terms <- df2$terms
ALL <- df2$Pr...t..
final2 <- data.frame(Terms, ALL)


## Location
var2 <- c("Influent", "Primary Clarifier Influent", "Primary Clarifier Effluent", "Final Clarifier Effluent", "UV Treated")
for (i in 1:length(var2)) {
  LocTerms <- paste0("ATE - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "UV Treated")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","UV - ",df$terms)
var2 <- c( "Final Clarifier Effluent", "Downstream")
for (i in 1:length(var2)) {
  LocTerms <- paste0("UV - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "Influent")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","INF - ",df$terms)
var2 <- c( "Residential Sewage", "Hospital Sewage")
for (i in 1:length(var2)) {
  LocTerms <- paste0("INF - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "Downstream")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","DSA - ",df$terms)
var2 <- c( "Upstream")
for (i in 1:length(var2)) {
  LocTerms <- paste0("DSA - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "Residential Sewage")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","RES - ",df$terms)
var2 <- c( "Hospital Sewage")
for (i in 1:length(var2)) {
  LocTerms <- paste0("RES - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "Primary Clarifier Influent")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","PCI - ",df$terms)
var2 <- c( "Primary Clarifier Effluent")
for (i in 1:length(var2)) {
  LocTerms <- paste0("PCI - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}


## Temperature
df$terms = gsub("Temperature","BT - ",df$terms)
df2 <- df[ df$terms %like% "BT" ,]
Terms <- df2$terms
ALL <- df2$Pr...t..
final2 <- rbind(final2, c(Terms, ALL))


## Media
df$terms = gsub("Media","LB - ",df$terms)
df2 <- df[ df$terms %like% "LB" ,]
Terms <- df2$terms
ALL <- df2$Pr...t..
final2 <- rbind(final2, c(Terms, ALL))

final <- merge(final2, final, by = "Terms")




## Antibiotic
frame$Antibiotic <- relevel(frame$Antibiotic, "Negative")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","HPC - ",df$terms)

var2 <- c( "Ampicillin", "Ciprofloxacin", "Doxycycline", "Sulfamethoxazole")
for (i in 1:length(var2)) {
  LocTerms <- paste0("HPC - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}


frame$Antibiotic <- relevel(frame$Antibiotic, "Ampicillin")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","AMP - ",df$terms)

var2 <- c( "Ciprofloxacin", "Doxycycline", "Sulfamethoxazole")
for (i in 1:length(var2)) {
  LocTerms <- paste0("AMP - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}



frame$Antibiotic <- relevel(frame$Antibiotic, "Ciprofloxacin")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","CIP - ",df$terms)

var2 <- c( "Doxycycline", "Sulfamethoxazole")
for (i in 1:length(var2)) {
  LocTerms <- paste0("CIP - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}



frame$Antibiotic <- relevel(frame$Antibiotic, "Doxycycline")
my.lm <- lm(logCFU_mL ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","DOX - ",df$terms)

var2 <- c( "Sulfamethoxazole")
for (i in 1:length(var2)) {
  LocTerms <- paste0("DOX - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}

order <- c("Mallard - Sugar", "BT - RT", "LB - R2A", "DSA - Upstream", "RES - Hospital Sewage", "INF - Residential Sewage", "INF - Hospital Sewage", "ATE - Influent", "PCI - Primary Clarifier Effluent", "ATE - Primary Clarifier Influent", "ATE - Primary Clarifier Effluent", "ATE - Final Clarifier Effluent", "ATE - UV Treated", "UV - Final Clarifier Effluent", "UV - Downstream", "HPC - Ampicillin", "HPC - Ciprofloxacin", "HPC - Doxycycline", "HPC - Sulfamethoxazole", "AMP - Ciprofloxacin", "AMP - Doxycycline", "AMP - Sulfamethoxazole", "CIP - Doxycycline", "CIP - Sulfamethoxazole", "DOX - Sulfamethoxazole")


table <- final %>%
  slice(match(order, Terms))

table$Terms = gsub("Upstream","UPA",table$Terms)
table$Terms = gsub("Hospital Sewage","HOS",table$Terms)
table$Terms = gsub("Residential Sewage","RES",table$Terms)
table$Terms = gsub("Primary Clarifier Influent","PCI",table$Terms)
table$Terms = gsub("Influent","INF",table$Terms)
table$Terms = gsub("Primary Clarifier Effluent","PCE",table$Terms)
table$Terms = gsub("Final Clarifier Effluent","FCE",table$Terms)
table$Terms = gsub("Downstream","DSA",table$Terms)
table$Terms = gsub("UV Treated","UV",table$Terms)
table$Terms = gsub("Ampicillin","AMP",table$Terms)
table$Terms = gsub("Ciprofloxacin","CIP",table$Terms)
table$Terms = gsub("Doxycycline","DOX",table$Terms)
table$Terms = gsub("Sulfamethoxazole","SUL",table$Terms)

table$Terms = gsub("DSA - UPA","UPA - DSA",table$Terms)
table$Terms = gsub("INF - RES","RES - INF",table$Terms)
table$Terms = gsub("INF - HOS","HOS - INF",table$Terms)
table$Terms = gsub("ATE - INF","INF - ATE",table$Terms)
table$Terms = gsub("ATE - PCI","PCI - ATE",table$Terms)
table$Terms = gsub("ATE - PCE","PCE - ATE",table$Terms)
table$Terms = gsub("UV - FCE","FCE - UV",table$Terms)

write.table(table, file=paste0(output, "Count_pValTable.tsv"), sep="\t",row.names=FALSE)

