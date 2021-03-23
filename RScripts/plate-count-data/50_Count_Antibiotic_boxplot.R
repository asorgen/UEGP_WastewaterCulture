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
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")

frame=read.table(paste0(input, "counts.txt"),sep="\t",header = TRUE)
frame=frame[frame$TimePoint=="4",]
frame$CFU <- frame$CFU.mL + 1

Antibiotic.aov <- aov(frame$CFU ~ frame$Antibiotic)
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

df$V1 = gsub("Ampicillin","1",df$V1)
df$V2 = gsub("Ampicillin","1",df$V2)

df$V1 = gsub("Ciprofloxacin","2",df$V1)
df$V2 = gsub("Ciprofloxacin","2",df$V2)

df$V1 = gsub("Doxycycline","3",df$V1)
df$V2 = gsub("Doxycycline","3",df$V2)

df$V1 = gsub("Sulfamethoxazole","4",df$V1)
df$V2 = gsub("Sulfamethoxazole","4",df$V2)

df$V1 = gsub("Negative","5",df$V1)
df$V2 = gsub("Negative","5",df$V2)

for (j in 1:length(df$sig.pvals)) {
  if (df$sig.pvals[j] < 0.001) {
    df[j,5] = "***"
  } else if (df$sig.pvals[j] > 0.001 & df$sig.pvals[j] < 0.01) {
    df[j,5] = "**"
  } else {
    df[j,5] = "*"
  }
}

capture.output(smry, df,file=paste0(output, "Antibiotic_count_ANOVA.txt"))

palette <- c("#5599E7","#9B35E8","#DF6B7B","#62B460","#EA9178")
antibiotic_labs=c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","No Antibiotic")
names(antibiotic_labs) <- c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","Negative")
frame$Antibiotic = factor(frame$Antibiotic, levels = c("Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole","Negative"))

plot <- ggplot(frame, aes(y=CFU, x=Antibiotic, fill=Antibiotic)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
  labs(x="Antibiotic", y=expression(paste("Abundance","\n",log[10],"(CFU/mL)")))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  geom_signif(annotation=formatC("***", digits=1), y_position=9.5, xmin=1, xmax=5, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("***", digits=1), y_position=9, xmin=2, xmax=5, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("***", digits=1), y_position=8.5, xmin=3, xmax=5, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("***", digits=1), y_position=8, xmin=4, xmax=5, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("*", digits=1), y_position=7, xmin=1, xmax=3, tip_length = c(0, 0),textsize = 6)+
  geom_signif(annotation=formatC("*", digits=1), y_position=6.5, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 6)+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  scale_y_log10()+
  stat_summary(fun=mean, geom="point", shape=16, size=5, color = "white")+
  theme(legend.position = "none")+
  scale_x_discrete(labels = antibiotic_labs)

ggsave(plot, file=paste0(output, "Count_Antibiotic_boxplot.pdf"),width = 20,height = 15)



