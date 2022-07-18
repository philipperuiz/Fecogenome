library(data.table)
library(tidyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(dplyr)
library(Hmisc)
library(plyr)
library(compositions)
library(zCompositions)
library(vegan)
library(scales)
library(viridisLite)
library(pals)
library(pheatmap)
library(pathview)

#all files can be downloaded at: https://doi.org/10.15454/VWZBUF

setwd("/home/panos/tmp/")
sample_info<-read.table("categories_samples2", sep="\t", header=TRUE, row.names = 1,quote="", fill=FALSE)
individuality=subset(sample_info, sample.type=="individual")
individuality=droplevels(individuality)




setwd("/home/panos/tmp/mapstat")
isos=read.table("all.genomic.corrected_n_imputed.fragment.by.length.ratios_130422")
isos5=read.table("all.genomic_NOT_imputed.fragment.by.length.ratios_130422")

isos2<-isos[!(row.names(isos) %in% "193782"), ]
listnames=rownames(individuality)
isos3 <- isos2 %>% dplyr::select(all_of(listnames))
isos6<-isos5[!(row.names(isos5) %in% "193782"), ]

isos7 <- isos6 %>% dplyr::select(all_of(listnames))
isos7=t(isos7)
isos7=as.data.frame(isos7)
filter <- isos7[,colSums(isos7) > 0]
names<-colnames(filter)
otu.table.filter <- subset(isos3, rownames(isos3) %in% names)
otu.table.filter<-droplevels(otu.table.filter)

#find most abundant and frequent contig
stats<-as.data.frame(rowSums(otu.table.filter))
temp<-as.data.frame(t(isos7))
sum(temp[2,]>0)
temp2<-subset(temp, rownames(temp) %in% rownames(stats))
temp3<-as.data.frame(matrix(ncol=1,nrow=nrow(temp2)))
rownames(temp3)<-rownames(temp2)
for (i in 1:nrow(temp2)) {
  temp3[i,]<-sum(temp2[i,]>0)}
stats2<-merge(stats,temp3,by="row.names") #temp3 is frequency, stats is abundnance
colnames(stats2)[2]<-c("totalcounts")
top<-stats2 %>% slice_max(totalcounts, n=1000)       # top 1000 in abundance 
top<-subset(top, V1 > 167)                                  # at least in 80% of samples cause it is 167/210
rownames(top)<-top$Row.names
extract2<-rownames(top)
extract2<-gsub("flag=.*.", "", extract2)
extract2<-gsub(".txt", ".txt_", extract2)


setwd("/home/panos/pCloudDrive/DELL-work/tmp-synced/")
setwd("/home/panos/tmp")
eggNOG <- read.table("MAGs.merged2.all.emapnotations.tbl", sep="\t", header=TRUE, row.names = 1, quote="", fill=FALSE)
eggNOG$names=rownames(eggNOG)
eggNOG$names<-sub("_[^_]+$", "", eggNOG$names)
eggNOG$mag<-eggNOG$names
eggNOG$mag<-gsub(".txt_.*","", eggNOG$mag)
eggNOG$mag<-gsub("(^.*)_W[0-9].*", "\\1", eggNOG$mag)
eggNOG$mag<-gsub("(^.*)_N[0-9].*", "\\1", eggNOG$mag)
data<-read.table("phylogeny_metainfo4v2",header=TRUE, row.names = 1, sep = ",")
data<-data[1:3]
colnames(data)[1]<-c("mag")
test7<-merge(eggNOG,data,by="mag",all.x = TRUE)
test8<-subset(test7, !(KEGG_ko == "-"))

##############################################################
eggNOG.selected.all<-test8 # # # # # # THIS IS ALL contigs 
eggNOG.selected.core<-subset(test8, names %in% extract2) # this is the contigs of the "core" MAGs
##############################################################




