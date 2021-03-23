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
frame$Condition <- paste(frame$Temperature, frame$Media, sep = "_")

Condition.aov <- aov(frame$shannon ~ frame$Condition)
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

df$V1 = gsub("BT_LB","2",df$V1)
df$V2 = gsub("BT_LB","2",df$V2)

df$V1 = gsub("RT_LB","3",df$V1)
df$V2 = gsub("RT_LB","3",df$V2)

df$V1 = gsub("BT_R2A","4",df$V1)
df$V2 = gsub("BT_R2A","4",df$V2)

df$V1 = gsub("RT_R2A","5",df$V1)
df$V2 = gsub("RT_R2A","5",df$V2)

df$V1 = gsub("_","1",df$V1)
df$V2 = gsub("_","1",df$V2)


frame$Location=factor(frame$Location)
frame$Condition=factor(frame$Condition)
frame$Temperature=factor(frame$Temperature)
frame$Media=factor(frame$Media)
frame$Site=factor(frame$Site)
palette <- c("#5599E7","#9B35E8","#DF6B7B","#62B460","#EA9178")


frame$Location=factor(frame$Location)
frame$Condition=factor(frame$Condition)
frame$Temperature=factor(frame$Temperature)
frame$Media=factor(frame$Media)
frame$Site=factor(frame$Site)

frame$Condition <- factor(frame$Condition, levels = c("n/a_n/a", "BT_LB", "RT_LB", "BT_R2A", "RT_R2A"))

df$stars <- "x"
for (j in 1:length(df$sig.pvals)) {
  if (df$sig.pvals[j] < 0.001) {
    df[j,5] = "***"
  } else if (df$sig.pvals[j] > 0.001 & df$sig.pvals[j] < 0.01) {
    df[j,5] = "**"
  } else {
    df[j,5] = "*"
  }
}

df

palette <- c("#5599E7","#9B35E8","#DF6B7B","#62B460","#EA9178")

plot <- ggplot(frame, aes(y=shannon, x=Condition, fill=Condition)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+
  stat_summary(fun=mean, colour="white", geom="point", shape=16, size=5)+
  labs(x="Condition", y="Shannon Diversity Index")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  geom_signif(annotation=formatC("***", digits=1), y_position=4.4, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 8)+
  geom_signif(annotation=formatC("***", digits=1), y_position=5.0, xmin=1, xmax=4, tip_length = c(0, 0),textsize = 8)+
  geom_signif(annotation=formatC("***", digits=1), y_position=4.7, xmin=1, xmax=3, tip_length = c(0, 0),textsize = 8)+
  geom_signif(annotation=formatC("***", digits=1), y_position=5.3, xmin=1, xmax=5, tip_length = c(0, 0),textsize = 8)+
  geom_signif(annotation=formatC("**", digits=1), y_position=4.0, xmin=3, xmax=4, tip_length = c(0, 0),textsize = 8)+
  geom_signif(annotation=formatC("**", digits=1), y_position=4.3, xmin=3, xmax=5, tip_length = c(0, 0),textsize = 8)+
  theme(legend.position = "none")+
  scale_x_discrete(labels = c("Environmental","LB agar at \n Body temp","LB agar at \n Room temp","R2A agar at \n Body temp","R2A agar at \n Room temp"))+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))

ggsave(plot, file=paste0(output, "Condition_SDI_boxplot.pdf"),width = 15,height = 10)




SampleType.t <- t.test(frame$shannon ~ frame$SampleType)

if (SampleType.t$p.value < 0.001) {
  star = "***"
} else if (SampleType.t$p.value > 0.001 & SampleType.t$p.value < 0.01) {
  star = "**"
} else {
  star = "*"
}

plot <- ggplot(frame, aes(y=shannon, x=SampleType, fill=SampleType)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+
  stat_summary(fun=mean, colour="white", geom="point", shape=16, size=5)+
  labs(x="SampleType", y="Shannon Diversity Index")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  geom_signif(annotation=star, y_position=4.4, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 6)+
  theme(legend.position = "none")+
  scale_x_discrete(labels = c("Culture","Environmental"))+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 24))+
  theme(axis.text.y = element_text(size = 24))

ggsave(plot, file=paste0(output, "SampleType_SDI_boxplot.pdf"),width = 20,height = 15)
