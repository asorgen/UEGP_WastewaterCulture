#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-18-21
#Description:

setwd('~/BioLockJ_pipelines/otu-analysis_2021Feb22/64_CountGraph/script/')
## Libraries
library(ggplot2)
library(grid)
library(scales)
library(ggsignif)
library(plyr)
library(plotrix)

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")


myT <- read.table(paste0(input, "counts.txt"), sep = "\t", header = TRUE)
T4 <- myT[myT$TimePoint=="4",]

T4$log.CFU.mL <- log10(T4$CFU.mL+1)
T4$CFU.mL <- T4$CFU.mL+1

data_summary <- function(data, varname, groupnames){
  require(plyr)
  require(plotrix)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm = TRUE),
      st.er = std.error(x[[col]], na.rm =TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun = summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df <- data_summary(T4, varname = "CFU.mL", groupnames = c("Location", "Antibiotic"))
# write.table(df,file = "Location_Antibiotic_count_summary.txt",row.names = FALSE,sep = "\t")
# df <- data_summary(T4, varname = "log.CFU.mL", groupnames = c("Location", "Antibiotic"))

Antibiotics <- factor(df$Antibiotic, levels = c("Negative","Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole"))
# df$Location=factor(df$Location, levels = c("Upstream","Residential Sewage","Hospital Sewage", "Influent", "Primary Clarifier Influent", "Primary Clarifier Effluent","Aeration Tank Effluent","Final Clarifier Effluent","UV Treated","Downstream"))
# Treatment_Stage=factor(df$Treatment, levels = c("Upstream","Pre-treatment","Post-treatment","Downstream"))


df$Location=factor(df$Location, levels = c("Upstream","Residential Sewage","Hospital Sewage", "Influent", "Primary Clarifier Influent", "Primary Clarifier Effluent","Aeration Tank Effluent","Final Clarifier Effluent","UV Treated","Downstream"))
df$Antibiotic=gsub("Negative","No Antibiotic",df$Antibiotic)
df$Antibiotic <- factor(df$Antibiotic, levels = c("No Antibiotic","Ampicillin","Ciprofloxacin","Doxycycline","Sulfamethoxazole"))
palette <- c("#5599E7","#9B35E8","#DF6B7B","#62B460","#EA9178")

plot <- ggplot(data = df, aes(x=Location, y=CFU.mL, colour = Antibiotic, group = Antibiotic))+
  geom_line(size=3)+
  geom_point(size=7)+
  labs(x="Treatment Stage", y=expression(paste("Abundance ",log[10],"(CFU/mL)")))+
  theme_bw()+
  theme(strip.text = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 35)) +
  theme(axis.title.x = element_text(size = 35))+
  theme(axis.text.x = element_text(size = 30))+
  theme(axis.text.y = element_text(size = 30))+
  theme(axis.title.x=element_text(margin=margin(t = 25, r = 0, b = 0, l = 0)))+
  geom_errorbar(aes(ymin=CFU.mL-st.er, ymax=CFU.mL+st.er), width=0.75, position = position_dodge(0))+
  theme(legend.title = element_blank(), legend.text = element_text(size = 30))+
  scale_color_manual(values = palette)+
  scale_x_discrete(labels = c("UP","RES","HOS","INF","PCI","PCE","ATE","FCE","UV","DS"))+
  #theme(legend.justification=c(1,0), legend.position=c(0.99,0.82))+
  theme(legend.box.background = element_rect(color="black", size=1))+
  scale_y_log10(breaks = c(1,1e1,1e2,1e3,1e4,1e5,1e6))

ggsave(plot, file=paste0(output, "Location_Antibiotic_CFU_graph.pdf"),width = 20,height = 10)


