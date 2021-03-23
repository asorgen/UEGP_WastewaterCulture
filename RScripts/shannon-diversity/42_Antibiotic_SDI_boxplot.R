#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: 

## Libraries
library(ggplot2)
library(grid)
library(scales)
library(ggsignif)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
inputModule = dir(pipeRoot, pattern="MetaUpdate", full.names=TRUE)
input = file.path(inputModule,"output/")
output = file.path(moduleDir,"output/")

frame=read.table(paste0(input, "metaUpdate.tsv"),sep="\t",header = TRUE)
frame <- frame[frame$SampleType == "Culture",]

Antibiotic.aov <- aov(frame$shannon ~ frame$Antibiotic)
smry <- summary.aov(Antibiotic.aov)

if (smry[[1]]$`Pr(>F)`[1] < 0.05) {
  tukey_Antibiotic <- TukeyHSD(x=Antibiotic.aov, conf.level=0.95)
  pvals <- tukey_Antibiotic$`frame$Antibiotic`[,4]
  
  sig.diffs <- vector()
  sig.pvals <- vector()
  index <- 1
  
  for (i in 1:length(pvals)) {
    if (pvals[i] < 0.05) {
      sig.diffs[index] <- names(pvals[i])
      sig.pvals[index] <- pvals[i]
      index <- index + 1
    }
  }
  df <- data.frame(sig.diffs, sig.pvals)
}

string <- strsplit(as.character(df$sig.diffs),split = "-")
temp_string=do.call(rbind,string)
df <- cbind(temp_string,df)
colnames(df)[colnames(df)=="1"] <- "V1"
colnames(df)[colnames(df)=="2"] <- "V2"

df$V1 = gsub("Amp","1",df$V1)
df$V2 = gsub("Amp","1",df$V2)

df$V1 = gsub("Cip","2",df$V1)
df$V2 = gsub("Cip","2",df$V2)

df$V1 = gsub("Dox","3",df$V1)
df$V2 = gsub("Dox","3",df$V2)

df$V1 = gsub("Sulf","4",df$V1)
df$V2 = gsub("Sulf","4",df$V2)

df$V1 = gsub("Neg","5",df$V1)
df$V2 = gsub("Neg","5",df$V2)


frame$Location=factor(frame$Location)
frame$Antibiotic=factor(frame$Antibiotic)
frame$Temperature=factor(frame$Temperature)
frame$Media=factor(frame$Media)
frame$Site=factor(frame$Site)
abx <- c("Neg","Amp","Cip","Dox","Sulf")
palette <- c("#5599E7","#9B35E8","#DF6B7B","#62B460","#EA9178")


frame$Location=factor(frame$Location)
frame$Antibiotic=factor(frame$Antibiotic)
frame$Temperature=factor(frame$Temperature)
frame$Media=factor(frame$Media)
frame$Site=factor(frame$Site)

frame$Antibiotic <- factor(frame$Antibiotic, levels = c("Amp","Cip","Dox","Sulf","Neg"))


for (j in 1:length(df$sig.pvals)) {
  if (df$sig.pvals[j] < 0.001) {
    df[j,5] = "***"
  } else if (df$sig.pvals[j] > 0.001 & df$sig.pvals[j] < 0.01) {
    df[j,5] = "**"
  } else {
    df[j,5] = "*"
  }
}
palette <- c("#5599E7","#9B35E8","#DF6B7B","#62B460","#EA9178")

plot <- ggplot(frame, aes(y=shannon, x=Antibiotic, fill=Antibiotic)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
  stat_summary(fun=mean, colour="white", geom="point", shape=16, size=5)+
  labs(x="Antibiotic", y="Shannon Diversity Index")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  geom_signif(annotation=formatC("***", digits=1), y_position=5.5, xmin=1, xmax=5, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("***", digits=1), y_position=5.2, xmin=2, xmax=5, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("*", digits=1), y_position=4.9, xmin=2, xmax=4, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("***", digits=1), y_position=4.6, xmin=3, xmax=5, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("**", digits=1), y_position=4.3, xmin=4, xmax=5, tip_length = c(0, 0),textsize = 6)+
  theme(legend.position = "none")+
  scale_x_discrete(labels = c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","No Antibiotic"))+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 24))+
  theme(axis.text.y = element_text(size = 24))
  
ggsave(plot, file=paste0(output, "Antibiotic_SDI_boxplot.pdf"),width = 20,height = 15)
