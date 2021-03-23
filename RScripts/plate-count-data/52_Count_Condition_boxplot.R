#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description: Generate colony count boxplot based on incubation conditions.

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
frame <- frame[frame$TimePoint=="4",]
frame$Location=factor(frame$Location)
frame$Antibiotic=factor(frame$Antibiotic, levels = c("Negative","Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole"))
frame$Temperature=factor(frame$Temperature)
frame$Media=factor(frame$Media)
frame$Site=factor(frame$Site)
frame$Condition <- paste(frame$Temperature,frame$Media,sep = "_")
frame$Condition <- factor(frame$Condition)
frame$CFU <- frame$CFU.mL + 1

Condition.aov <- aov(frame$CFU ~ frame$Condition)
smry <- summary.aov(Condition.aov)

if (smry[[1]]$`Pr(>F)`[1] < 0.05) {
  tukey_Condition <- TukeyHSD(x=Condition.aov, conf.level=0.95)
  pvals <- tukey_Condition$`frame$Condition`[,4]
  
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

df$V1 = gsub("BT_LB","1",df$V1)
df$V2 = gsub("BT_LB","1",df$V2)

df$V1 = gsub("RT_LB","2",df$V1)
df$V2 = gsub("RT_LB","2",df$V2)

df$V1 = gsub("BT_R2A","3",df$V1)
df$V2 = gsub("BT_R2A","3",df$V2)

df$V1 = gsub("RT_R2A","4",df$V1)
df$V2 = gsub("RT_R2A","4",df$V2)

for (j in 1:length(df$sig.pvals)) {
  if (df$sig.pvals[j] < 0.001) {
    df[j,5] = "***"
  } else if (df$sig.pvals[j] > 0.001 & df$sig.pvals[j] < 0.01) {
    df[j,5] = "**"
  } else {
    df[j,5] = "*"
  }
}

capture.output(smry, df,file=paste0(output, "Condition_count_ANOVA.txt"))

palette <- c("#5599E7","#9B35E8","#DF6B7B","#62B460","#EA9178")
antibiotic_labs=c("LB agar\nBody temp","LB agar\nRoom temp","R2A agar\nBody temp","R2A agar\nRoom temp")
names(antibiotic_labs) <- c("BT_LB","RT_LB","BT_R2A","RT_R2A")
frame$Condition = factor(frame$Condition, levels = c("BT_LB","RT_LB","BT_R2A","RT_R2A"))

plot <- ggplot(frame, aes(y=CFU, x=Condition, fill=Condition)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
  labs(x="Condition", y=expression(paste("Abundance","\n",log[10],"(CFU/mL)")))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  geom_signif(annotation=formatC(paste0(df[1,5]), digits=1), y_position=7.5, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 8)+
  geom_signif(annotation=formatC("***", digits=1), y_position=8, xmin=1, xmax=4, tip_length = c(0, 0),textsize = 8)+
  geom_signif(annotation=formatC("**", digits=1), y_position=7.5, xmin=3, xmax=4, tip_length = c(0, 0),textsize = 8)+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  scale_y_log10()+
  stat_summary(fun=mean, geom="point", shape=16, size=5, color = "white")+
  theme(legend.position = "none")+
  scale_x_discrete(labels = antibiotic_labs)

ggsave(plot, file=paste0(output, "Count_Condition_boxplot.pdf"),width = 15,height = 10)






Temperature.t <- t.test(frame$CFU ~ frame$Temperature)
capture.output(Temperature.t,file=paste0(output, "Temperature_count_tTest.txt"))
pVal <- Temperature.t$p.value

if (pVal < 0.001) {
  stars = "***"
} else if (pVal > 0.001 & pVal < 0.01) {
  stars = "**"
} else if (pVal > 0.01 & pVal < 0.05) {
  stars = "*"
} else {
  stars = "NS"
}

antibiotic_labs=c("Body temp","Room temp")
names(antibiotic_labs) <- c("BT","RT")

plot <- ggplot(frame, aes(y=CFU, x=Temperature, fill=Temperature)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
  labs(x="Temperature", y=expression(paste("Abundance","\n",log[10],"(CFU/mL)")))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  geom_signif(annotation=formatC(paste0(stars), digits=1), y_position=7.5, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 6)+
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.title.x = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  scale_y_log10()+
  stat_summary(fun=mean, geom="point", shape=16, size=5, color = "white")+
  theme(legend.position = "none")+
  scale_x_discrete(labels = antibiotic_labs)

ggsave(plot, file=paste0(output, "Count_Temperature_boxplot.pdf"),width = 15,height = 10)





Media.t <- t.test(frame$CFU ~ frame$Media)
capture.output(Media.t,file=paste0(output, "Media_count_tTest.txt"))
pVal <- Media.t$p.value

if (pVal < 0.001) {
  stars = "***"
} else if (pVal > 0.001 & pVal < 0.01) {
  stars = "**"
} else if (pVal > 0.01 & pVal < 0.05) {
  stars = "*"
} else {
  stars = "NS"
}


antibiotic_labs=c("LB agar","R2A agar")
names(antibiotic_labs) <- c("LB","R2A")

plot <- ggplot(frame, aes(y=CFU, x=Media, fill=Media)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
  labs(x="Media", y=expression(paste("Abundance","\n",log[10],"(CFU/mL)")))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  geom_signif(annotation=formatC(paste0(stars), digits=1), y_position=7.5, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 6)+
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.title.x = element_text(size = 30))+
  theme(axis.text.x = element_text(size = 25))+
  theme(axis.text.y = element_text(size = 25))+
  scale_y_log10()+
  stat_summary(fun=mean, geom="point", shape=16, size=5, color = "white")+
  theme(legend.position = "none")+
  scale_x_discrete(labels = antibiotic_labs)

ggsave(plot, file=paste0(output, "Count_Media_boxplot.pdf"),width = 15,height = 10)

