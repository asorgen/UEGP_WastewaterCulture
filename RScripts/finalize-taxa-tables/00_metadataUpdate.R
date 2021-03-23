#Author: Alicia Sorgen
#BioLockJ configuration: Alicia Sorgen
#Date: 01-05-21
#Description: Update metadata file to include Shannon diversity and observed OTU data

rm(list=ls())
pipeRoot = dirname(dirname(getwd()))
moduleDir <- dirname(getwd())
input = file.path(pipeRoot,"input/")
output = file.path(moduleDir,"output/")

metadata <- read.table(paste0(input, "metadata.tsv"), sep = "\t", header = TRUE)

shannon <- read.table(paste0(input, "shannon.tsv"), sep = "\t", header = TRUE)
names(shannon)[names(shannon)=="Sample.ID"]="SampleID"

observed_otus <- read.table(paste0(input, "observed_otus.tsv"), sep = "\t", header = TRUE)
names(observed_otus)[names(observed_otus)=="Sample.ID"]="SampleID"

seqs=read.table(paste0(input, "SEQfreqperSample.csv"),sep=",",header = FALSE)
colnames(seqs)[colnames(seqs)=="V1"] <- "SampleID"
colnames(seqs)[colnames(seqs)=="V2"] <- "SequenceFrequency"

merge=merge(metadata,shannon,by="SampleID")
merge=merge(merge,observed_otus,by="SampleID")
merge=merge(merge,seqs,by="SampleID")

write.table(x = merge, file = paste0(output, "metaUpdate.tsv"), row.names = FALSE, sep="\t")
