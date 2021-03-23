#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: 

## Libraries
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(plyr)

rm(list=ls())
##### Prep #####
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputnModule = dir(pipeRoot, pattern="MetaUpdate", full.names=TRUE)
inputPath = file.path(inputnModule,"output/")
output = file.path(moduleDir,"output/")

frame=read.table(paste0(inputPath, "metaUpdate.tsv"),sep="\t",header = TRUE)
frame <- frame[frame$SampleType == "Culture",]

frame$Site <- factor(frame$Site)
frame$Location <- factor(frame$Location)
frame$Temperature <- factor(frame$Temperature)
frame$Media <- factor(frame$Media)
frame$Antibiotic <- factor(frame$Antibiotic)

##### ARB p values #####
frame2 <- frame[!(frame$Antibiotic == "Neg"),]
frame2$Site <- factor(frame2$Site)
frame2$Location <- factor(frame2$Location)
frame2$Temperature <- factor(frame2$Temperature)
frame2$Media <- factor(frame2$Media)
frame2$Antibiotic <- factor(frame2$Antibiotic)


my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
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
var2 <- c("INF", "PCI", "PCE", "FCE", "UV")
for (i in 1:length(var2)) {
  LocTerms <- paste0("ATE - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "UV")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","UV - ",df$terms)
var2 <- c( "FCE", "DSA")
for (i in 1:length(var2)) {
  LocTerms <- paste0("UV - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "INF")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","INF - ",df$terms)
var2 <- c( "RES", "HOS")
for (i in 1:length(var2)) {
  LocTerms <- paste0("INF - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "DSA")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","DSA - ",df$terms)
var2 <- c( "UPA")
for (i in 1:length(var2)) {
  LocTerms <- paste0("DSA - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "RES")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","RES - ",df$terms)
var2 <- c( "HOS")
for (i in 1:length(var2)) {
  LocTerms <- paste0("RES - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ARB <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ARB))
  
}

frame2$Location <- relevel(frame2$Location, "PCI")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame2)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","PCI - ",df$terms)
var2 <- c( "PCE")
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


##### Individual antibiotic p values #####
var1 <- c("Neg", "Amp", "Cip", "Dox", "Sulf")

for (j in 1:length(var1)) {
  Antibiotic <- paste0(var1[j])
  frame2 <- frame[frame$Antibiotic %in% Antibiotic,]
  frame2$Site <- factor(frame2$Site)
  frame2$Location <- factor(frame2$Location)
  frame2$Temperature <- factor(frame2$Temperature)
  frame2$Media <- factor(frame2$Media)
  
  
  my.lm <- lm(shannon ~ Site + Location + Temperature + Media, data = frame2)
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
  var2 <- c("INF", "PCI", "PCE", "FCE", "UV")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("ATE - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  counts <- ldply(frame2$Location, function(c) sum(c=="UV"))  
  
  if (sum(counts$V1) < 1) {
    
    var2 <- c( "FCE", "DSA")
    for (i in 1:length(var2)) {
      Terms <- paste0("UV - ", var2[i])
      ARB <- 0
      final2 <- rbind(final2, c(Terms, ARB))
    }
    
  } else {
    frame2$Location <- relevel(frame2$Location, "UV")
    my.lm <- lm(shannon ~ Site + Location + Temperature + Media , data = frame2)
    smry <- summary(my.lm)
    coef <- smry$coefficients
    terms <- rownames(coef)
    df <-cbind(terms, coef)
    df <- data.frame(df)
    df$terms = gsub("Location","UV - ",df$terms)
    var2 <- c( "FCE", "DSA")
    for (i in 1:length(var2)) {
      LocTerms <- paste0("UV - ", var2[i])
      df2 <- df[ df$terms %in% LocTerms ,]
      Terms <- df2$terms
      ARB <- df2$Pr...t..
      final2 <- rbind(final2, c(Terms, ARB))
    }
    
  }

  
  
  frame2$Location <- relevel(frame2$Location, "INF")
  my.lm <- lm(shannon ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","INF - ",df$terms)
  var2 <- c( "RES", "HOS")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("INF - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  counts <- ldply(frame2$Location, function(c) sum(c=="DSA"))  
  
  if (sum(counts$V1) < 1) {
    
    var2 <- c( "UPA")
    for (i in 1:length(var2)) {
      Terms <- paste0("DSA - ", var2[i])
      ARB <- 0
      final2 <- rbind(final2, c(Terms, ARB))
    }
    
  } else {
    
    frame2$Location <- relevel(frame2$Location, "DSA")
    my.lm <- lm(shannon ~ Site + Location + Temperature + Media , data = frame2)
    smry <- summary(my.lm)
    coef <- smry$coefficients
    terms <- rownames(coef)
    df <-cbind(terms, coef)
    df <- data.frame(df)
    df$terms = gsub("Location","DSA - ",df$terms)
    var2 <- c( "UPA")
    for (i in 1:length(var2)) {
      LocTerms <- paste0("DSA - ", var2[i])
      df2 <- df[ df$terms %in% LocTerms ,]
      Terms <- df2$terms
      ARB <- df2$Pr...t..
      final2 <- rbind(final2, c(Terms, ARB))
    }
    
  }
  
  frame2$Location <- relevel(frame2$Location, "RES")
  my.lm <- lm(shannon ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","RES - ",df$terms)
  var2 <- c( "HOS")
  for (i in 1:length(var2)) {
    LocTerms <- paste0("RES - ", var2[i])
    df2 <- df[ df$terms %in% LocTerms ,]
    Terms <- df2$terms
    ARB <- df2$Pr...t..
    final2 <- rbind(final2, c(Terms, ARB))
    
  }
  
  frame2$Location <- relevel(frame2$Location, "PCI")
  my.lm <- lm(shannon ~ Site + Location + Temperature + Media , data = frame2)
  smry <- summary(my.lm)
  coef <- smry$coefficients
  terms <- rownames(coef)
  df <-cbind(terms, coef)
  df <- data.frame(df)
  df$terms = gsub("Location","PCI - ",df$terms)
  var2 <- c( "PCE")
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
  colnames(final2)[colnames(final2)=="ARB"] <- Antibiotic
  
  final <- merge(final, final2, by = "Terms")
  
}



##### ALL p values #####
frame$Site <- factor(frame$Site)
frame$Location <- factor(frame$Location)
frame$Temperature <- factor(frame$Temperature)
frame$Media <- factor(frame$Media)
frame$Antibiotic <- factor(frame$Antibiotic)


my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
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
var2 <- c("INF", "PCI", "PCE", "FCE", "UV")
for (i in 1:length(var2)) {
  LocTerms <- paste0("ATE - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "UV")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","UV - ",df$terms)
var2 <- c( "FCE", "DSA")
for (i in 1:length(var2)) {
  LocTerms <- paste0("UV - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "INF")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","INF - ",df$terms)
var2 <- c( "RES", "HOS")
for (i in 1:length(var2)) {
  LocTerms <- paste0("INF - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "DSA")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","DSA - ",df$terms)
var2 <- c( "UPA")
for (i in 1:length(var2)) {
  LocTerms <- paste0("DSA - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "RES")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","RES - ",df$terms)
var2 <- c( "HOS")
for (i in 1:length(var2)) {
  LocTerms <- paste0("RES - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  final2 <- rbind(final2, c(Terms, ALL))
  
}

frame$Location <- relevel(frame$Location, "PCI")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Location","PCI - ",df$terms)
var2 <- c( "PCE")
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
frame$Antibiotic <- relevel(frame$Antibiotic, "Neg")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","HPC - ",df$terms)

var2 <- c( "Amp", "Cip", "Dox", "Sulf")
for (i in 1:length(var2)) {
  LocTerms <- paste0("HPC - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}


frame$Antibiotic <- relevel(frame$Antibiotic, "Amp")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","AMP - ",df$terms)

var2 <- c( "Cip", "Dox", "Sulf")
for (i in 1:length(var2)) {
  LocTerms <- paste0("AMP - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}



frame$Antibiotic <- relevel(frame$Antibiotic, "Cip")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","CIP - ",df$terms)

var2 <- c( "Dox", "Sulf")
for (i in 1:length(var2)) {
  LocTerms <- paste0("CIP - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}



frame$Antibiotic <- relevel(frame$Antibiotic, "Dox")
my.lm <- lm(shannon ~ Site + Location + Temperature + Media + Antibiotic, data = frame)
smry <- summary(my.lm)
coef <- smry$coefficients
terms <- rownames(coef)
df <-cbind(terms, coef)
df <- data.frame(df)
df$terms = gsub("Antibiotic","DOX - ",df$terms)

var2 <- c( "Sulf")
for (i in 1:length(var2)) {
  LocTerms <- paste0("DOX - ", var2[i])
  df2 <- df[ df$terms %in% LocTerms ,]
  Terms <- df2$terms
  ALL <- df2$Pr...t..
  new <- c(Terms, ALL, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a")
  final <- rbind(final, new)
  
}

order <- c("Mallard - Sugar", "BT - RT", "LB - R2A", "DSA - UPA", "RES - HOS", "INF - RES", "INF - HOS", "ATE - INF", "PCI - PCE", "ATE - PCI", "ATE - PCE", "ATE - FCE", "ATE - UV", "UV - FCE", "UV - DSA", "HPC - Amp", "HPC - Cip", "HPC - Dox", "HPC - Sulf", "AMP - Cip", "AMP - Dox", "AMP - Sulf", "CIP - Dox", "CIP - Sulf", "DOX - Sulf")


table <- final %>%
  slice(match(order, Terms))


table$Terms = gsub("Amp","AMP",table$Terms)
table$Terms = gsub("Cip","CIP",table$Terms)
table$Terms = gsub("Dox","DOX",table$Terms)
table$Terms = gsub("Sulf","SUL",table$Terms)

table$Terms = gsub("DSA - UPA","UPA - DSA",table$Terms)
table$Terms = gsub("INF - RES","RES - INF",table$Terms)
table$Terms = gsub("INF - HOS","HOS - INF",table$Terms)
table$Terms = gsub("ATE - INF","INF - ATE",table$Terms)
table$Terms = gsub("ATE - PCI","PCI - ATE",table$Terms)
table$Terms = gsub("ATE - PCE","PCE - ATE",table$Terms)
table$Terms = gsub("UV - FCE","FCE - UV",table$Terms)

write.table(table, file=paste0(output, "SDI_pValTable.tsv"), sep="\t",row.names=FALSE)
