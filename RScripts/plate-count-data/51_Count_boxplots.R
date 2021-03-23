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

myT=read.table(paste0(input, "counts.txt"),sep="\t",header = TRUE)
myT=myT[myT$TimePoint=="4",]


df=myT[myT$Location=="Upstream"|myT$Location=="Downstream",]
df$Site=factor(df$Site)
df$Location=factor(df$Location, levels = c("Upstream","Downstream"))
df$Antibiotic=factor(df$Antibiotic, levels = c("Negative","Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole"))
df$Temperature=factor(df$Temperature)
df$Media=factor(df$Media)

myLm_ALL <- lm(formula = log10(CFU.mL+1) ~ Location, data = df)
summary(myLm_ALL)
myAnova_ALL <- anova(myLm_ALL)
myAnova_ALL

HPC <- df[df$Antibiotic == "Negative",]
ARB <- df[!(df$Antibiotic == "Negative"),]

myLm_HPC <- lm(formula = log10(CFU.mL+1) ~ Location, data = HPC)
myAnova_HPC <- anova(myLm_HPC)
myAnova_HPC
myLm_ARB <- lm(formula = log10(CFU.mL+1) ~ Location, data = ARB)
myAnova_ARB <- anova(myLm_ARB)
myAnova_ARB

capture.output(myAnova_ALL,myAnova_HPC,myAnova_ARB,file=paste0(output, "UPAvDSA_count_GLM.txt"))


ARB$Ab_type <- "Combined ARB"
HPC$Ab_type <- "Total heterotrophic growth"
df <- rbind(ARB, HPC)
df$new=df$CFU.mL+1

plot <- ggplot(df, aes(y=log10(new), x=Location)) +
  geom_boxplot()+
  labs(x="", y="Abundance Log10(CFU/mL)")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_minimal()+
  stat_summary(fun=mean, color="red", geom="point", shape=18, position=position_dodge(.73), size=5)+
  geom_signif(annotation = c("NS"), y_position = c(7.2), xmin = c(1), xmax = c(2), tip_length = c(0, 0),textsize = 6)+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 20))+
  theme(axis.text.y = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))+
  facet_grid(cols = vars(Ab_type))+
  theme(strip.text.x = element_text(size = 20, colour = "black"))+
  theme(legend.position = "none")
ggsave(plot, file=paste0(output, "UPAvDSA_count_boxplot.pdf"), width = 10, height = 7)






df=myT[myT$Location=="Residential Sewage"|myT$Location=="Hospital Sewage",]
df$Site=factor(df$Site)
df$Location=factor(df$Location, levels = c("Residential Sewage","Hospital Sewage"))
df$Antibiotic=factor(df$Antibiotic, levels = c("Negative","Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole"))
df$Temperature=factor(df$Temperature)
df$Media=factor(df$Media)

myLm_ALL <- lm(formula = log10(CFU.mL+1) ~ Location, data = df)
myAnova_ALL <- anova(myLm_ALL)
myAnova_ALL

HPC <- df[df$Antibiotic == "Negative",]
ARB <- df[!(df$Antibiotic == "Negative"),]

myLm_HPC <- lm(formula = log10(CFU.mL+1) ~ Location, data = HPC)
myAnova_HPC <- anova(myLm_HPC)
myAnova_HPC
myLm_ARB <- lm(formula = log10(CFU.mL+1) ~ Location, data = ARB)
myAnova_ARB <- anova(myLm_ARB)
myAnova_ARB

capture.output(myAnova_ALL,myAnova_HPC,myAnova_ARB,file=paste0(output, "RESvHOS_count_GLM.txt"))

ARB$Ab_type <- "Combined ARB"
HPC$Ab_type <- "Total heterotrophic growth"
df <- rbind(ARB, HPC)
df$new=df$CFU.mL+1

plot <- ggplot(df, aes(y=log10(new), x=Location)) +
  geom_boxplot()+
  labs(x="", y="Abundance Log10(CFU/mL)")+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_minimal()+
  stat_summary(fun=mean, color="red", geom="point", shape=18, position=position_dodge(.73), size=5)+
  geom_signif(annotation = c("*"), y_position = c(7.2), xmin = c(1), xmax = c(2), tip_length = c(0, 0),textsize = 6)+
  # scale_fill_brewer(palette="Paired")+
  theme(axis.title.y = element_text(size = 25)) +
  theme(axis.title.x = element_text(size = 25))+
  theme(axis.text.x = element_text(size = 18))+
  theme(axis.text.y = element_text(size = 20))+
  facet_grid(cols = vars(Ab_type))+
  theme(strip.text.x = element_text(size = 20, colour = "black"))+
  scale_x_discrete(labels = c("Residential\nSewage","Hospital\nSewage"))+
  theme(legend.position = "none")
ggsave(plot, file=paste0(output, "RESvHOS_count_boxplot.pdf"), width = 10, height = 7)




