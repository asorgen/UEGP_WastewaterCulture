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
frame <- frame[!(frame$Antibiotic == "Neg"),]

Location.aov <- aov(frame$shannon ~ frame$Location)
smry <- summary.aov(Location.aov)

# if (smry[[1]]$`Pr(>F)`[1] < 0.05) {
#   tukey_Location <- TukeyHSD(x=Location.aov, conf.level=0.95)
#   pvals <- tukey_Location$`frame$Location`[,4]
#   
#   sig.diffs <- vector()
#   sig.pvals <- vector()
#   index <- 1
#   
#   for (i in 1:length(pvals)) {
#     if (pvals[i] < 0.05) {
#       sig.diffs[index] <- names(pvals[i])
#       sig.pvals[index] <- pvals[i]
#       index <- index + 1
#     }
#   }
# 
# }


frame$Location <- factor(frame$Location, levels = c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))

palette <- c("#5599E7","#9B35E8","#DF6B7B","#9CF270","#EA9178","#C7F088","#E25047","#CBA7F2","#9DCBCB","#EDCDD0")


plot <- ggplot(frame, aes(y=shannon, x=Location, fill=Location)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
  stat_summary(fun=mean, colour="white", geom="point", shape=16, size=5)+
  labs(x="Location", y="Shannon Diversity Index")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  scale_fill_manual(values=palette) + theme_minimal()+
  # geom_signif(annotation=formatC("***", digits=1), y_position=4.4, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 6)+
  # geom_signif(annotation=formatC("***", digits=1), y_position=5.0, xmin=1, xmax=4, tip_length = c(0, 0),textsize = 6)+
  # geom_signif(annotation=formatC("***", digits=1), y_position=4.7, xmin=1, xmax=3, tip_length = c(0, 0),textsize = 6)+
  # geom_signif(annotation=formatC("***", digits=1), y_position=5.3, xmin=1, xmax=5, tip_length = c(0, 0),textsize = 6)+
  # geom_signif(annotation=formatC("**", digits=1), y_position=4.0, xmin=3, xmax=4, tip_length = c(0, 0),textsize = 6)+
  # geom_signif(annotation=formatC("**", digits=1), y_position=4.3, xmin=3, xmax=5, tip_length = c(0, 0),textsize = 6)+
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 24))+
  theme(axis.text.y = element_text(size = 24))

ggsave(plot, file=paste0(output, "Location_SDI_boxplot_ARBonly.pdf"),width = 20,height = 15)








# frame=read.table(paste0(input, "metaUpdate.tsv"),sep="\t",header = TRUE)
# frame <- frame[frame$SampleType == "Culture",]
# frame <- frame[(frame$Antibiotic == "Neg"),]
# 
# Location.aov <- aov(frame$shannon ~ frame$Location)
# smry <- summary.aov(Location.aov)
# smry
# 
# if (smry[[1]]$`Pr(>F)`[1] < 0.05) {
#   tukey_Location <- TukeyHSD(x=Location.aov, conf.level=0.95)
#   pvals <- tukey_Location$`frame$Location`[,4]
#   
#   sig.diffs <- vector()
#   sig.pvals <- vector()
#   index <- 1
#   
#   for (i in 1:length(pvals)) {
#     if (pvals[i] < 0.05) {
#       sig.diffs[index] <- names(pvals[i])
#       sig.pvals[index] <- pvals[i]
#       index <- index + 1
#     }
#   }
#   df <- data.frame(sig.diffs, sig.pvals)
#   string <- strsplit(as.character(df$sig.diffs),split = "-")
#   temp_string=do.call(rbind,string)
#   df <- cbind(temp_string,df)
#   colnames(df)[colnames(df)=="1"] <- "V1"
#   colnames(df)[colnames(df)=="2"] <- "V2"
#   
#   df$V1 = gsub("UPA","1",df$V1)
#   df$V2 = gsub("UPA","1",df$V2)
#   
#   df$V1 = gsub("RES","2",df$V1)
#   df$V2 = gsub("RES","2",df$V2)
#   
#   df$V1 = gsub("HOS","3",df$V1)
#   df$V2 = gsub("HOS","3",df$V2)
#   
#   df$V1 = gsub("INF","4",df$V1)
#   df$V2 = gsub("INF","4",df$V2)
#   
#   df$V1 = gsub("PCI","5",df$V1)
#   df$V2 = gsub("PCI","5",df$V2)
#   
#   df$V1 = gsub("PCE","6",df$V1)
#   df$V2 = gsub("PCE","6",df$V2)
#   
#   df$V1 = gsub("ATE","7",df$V1)
#   df$V2 = gsub("ATE","7",df$V2)
#   
#   df$V1 = gsub("FCE","8",df$V1)
#   df$V2 = gsub("FCE","8",df$V2)
#   
#   df$V1 = gsub("UV","9",df$V1)
#   df$V2 = gsub("UV","9",df$V2)
#   
#   df$V1 = gsub("DSA","10",df$V1)
#   df$V2 = gsub("DSA","10",df$V2)
#   
#   df$stars <- "x"
#   for (j in 1:length(df$sig.pvals)) {
#     if (df$sig.pvals[j] < 0.001) {
#       df[j,5] = "***"
#     } else if (df$sig.pvals[j] > 0.001 & df$sig.pvals[j] < 0.01) {
#       df[j,5] = "**"
#     } else {
#       df[j,5] = "*"
#     }
#   }
#   
#   df
#   
#   
# }
# 
# 
# frame$Location <- factor(frame$Location, levels = c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))
# 
# palette <- c("#5599E7","#9B35E8","#DF6B7B","#9CF270","#EA9178","#C7F088","#E25047","#CBA7F2","#9DCBCB","#EDCDD0")
# 
# 
# plot <- ggplot(frame, aes(y=shannon, x=Location, fill=Location)) +
#   geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
#   stat_summary(fun=mean, colour="white", geom="point", shape=16, size=5)+
#   labs(x="Location", y="Shannon Diversity Index")+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   scale_fill_manual(values=palette) + theme_minimal()+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=4.4, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=5.0, xmin=1, xmax=4, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=4.7, xmin=1, xmax=3, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=5.3, xmin=1, xmax=5, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("**", digits=1), y_position=4.0, xmin=3, xmax=4, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("**", digits=1), y_position=4.3, xmin=3, xmax=5, tip_length = c(0, 0),textsize = 6)+
#   theme(legend.position = "none")+
#   theme(axis.title.y = element_text(size = 25)) +
#   theme(axis.title.x = element_text(size = 25))+
#   theme(axis.text.x = element_text(size = 24))+
#   theme(axis.text.y = element_text(size = 24))
# 
# ggsave(plot, file=paste0(output, "Location_SDI_boxplot_HPConly.pdf"),width = 20,height = 15)
# 
# 
# 
# 
# 
# 
# frame=read.table(paste0(input, "metaUpdate.tsv"),sep="\t",header = TRUE)
# frame <- frame[frame$SampleType == "Culture",]
# 
# Location.aov <- aov(frame$shannon ~ frame$Location)
# smry <- summary.aov(Location.aov)
# smry
# 
# if (smry[[1]]$`Pr(>F)`[1] < 0.05) {
#   tukey_Location <- TukeyHSD(x=Location.aov, conf.level=0.95)
#   pvals <- tukey_Location$`frame$Location`[,4]
#   
#   sig.diffs <- vector()
#   sig.pvals <- vector()
#   index <- 1
#   
#   for (i in 1:length(pvals)) {
#     if (pvals[i] < 0.05) {
#       sig.diffs[index] <- names(pvals[i])
#       sig.pvals[index] <- pvals[i]
#       index <- index + 1
#     }
#   }
#   df <- data.frame(sig.diffs, sig.pvals)
#   string <- strsplit(as.character(df$sig.diffs),split = "-")
#   temp_string=do.call(rbind,string)
#   df <- cbind(temp_string,df)
#   colnames(df)[colnames(df)=="1"] <- "V1"
#   colnames(df)[colnames(df)=="2"] <- "V2"
#   
#   df$V1 = gsub("UPA","1",df$V1)
#   df$V2 = gsub("UPA","1",df$V2)
#   
#   df$V1 = gsub("RES","2",df$V1)
#   df$V2 = gsub("RES","2",df$V2)
#   
#   df$V1 = gsub("HOS","3",df$V1)
#   df$V2 = gsub("HOS","3",df$V2)
#   
#   df$V1 = gsub("INF","4",df$V1)
#   df$V2 = gsub("INF","4",df$V2)
#   
#   df$V1 = gsub("PCI","5",df$V1)
#   df$V2 = gsub("PCI","5",df$V2)
#   
#   df$V1 = gsub("PCE","6",df$V1)
#   df$V2 = gsub("PCE","6",df$V2)
#   
#   df$V1 = gsub("ATE","7",df$V1)
#   df$V2 = gsub("ATE","7",df$V2)
#   
#   df$V1 = gsub("FCE","8",df$V1)
#   df$V2 = gsub("FCE","8",df$V2)
#   
#   df$V1 = gsub("UV","9",df$V1)
#   df$V2 = gsub("UV","9",df$V2)
#   
#   df$V1 = gsub("DSA","10",df$V1)
#   df$V2 = gsub("DSA","10",df$V2)
#   
#   df$stars <- "x"
#   for (j in 1:length(df$sig.pvals)) {
#     if (df$sig.pvals[j] < 0.001) {
#       df[j,5] = "***"
#     } else if (df$sig.pvals[j] > 0.001 & df$sig.pvals[j] < 0.01) {
#       df[j,5] = "**"
#     } else {
#       df[j,5] = "*"
#     }
#   }
#   
#   df
#   
#   
# }
# 
# 
# frame$Location <- factor(frame$Location, levels = c("UPA","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DSA"))
# 
# palette <- c("#5599E7","#9B35E8","#DF6B7B","#9CF270","#EA9178","#C7F088","#E25047","#CBA7F2","#9DCBCB","#EDCDD0")
# 
# 
# plot <- ggplot(frame, aes(y=shannon, x=Location, fill=Location)) +
#   geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=2)+ 
#   stat_summary(fun=mean, colour="white", geom="point", shape=16, size=5)+
#   labs(x="Location", y="Shannon Diversity Index")+
#   geom_jitter(shape=16, position=position_jitter(0.2))+
#   scale_fill_manual(values=palette) + theme_minimal()+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=4.4, xmin=1, xmax=2, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=5.0, xmin=1, xmax=4, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=4.7, xmin=1, xmax=3, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("***", digits=1), y_position=5.3, xmin=1, xmax=5, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("**", digits=1), y_position=4.0, xmin=3, xmax=4, tip_length = c(0, 0),textsize = 6)+
#   # geom_signif(annotation=formatC("**", digits=1), y_position=4.3, xmin=3, xmax=5, tip_length = c(0, 0),textsize = 6)+
#   theme(legend.position = "none")+
#   theme(axis.title.y = element_text(size = 25)) +
#   theme(axis.title.x = element_text(size = 25))+
#   theme(axis.text.x = element_text(size = 24))+
#   theme(axis.text.y = element_text(size = 24))
# 
# ggsave(plot, file=paste0(output, "Location_SDI_boxplot.pdf"),width = 20,height = 15)
# 
