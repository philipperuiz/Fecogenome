library(plyr)
library(reshape2)
library(dplyr)
library(tidyr)
library(ape)
library(vegan)
library(MASS)
library(data.table)
library(indicspecies)
library(tidyverse)
library(Hmisc)
library(compositions)
library(zCompositions)
library(pals)
library(pheatmap)
library(pathview)

#all files can be downloaded at: https://doi.org/10.15454/VWZBUF

setwd("/home/panos/tmp/")
sample_info<-read.table("categories_samples2", sep="\t", header=TRUE, row.names = 1,quote="", fill=FALSE)
individuality=subset(sample_info, sample.type=="individual")
individuality=droplevels(individuality)

individuality2<-individuality[complete.cases(individuality$diet),]
individuality=droplevels(individuality)
individuality2=droplevels(individuality2)

eggNOG <- read.table("MAGs.merged2.all.emapnotations.tbl", sep="\t", header=TRUE, row.names = 1, quote="", fill=FALSE)
eggNOG$names=rownames(eggNOG)
eggNOG$names<-sub("_[^_]+$", "", eggNOG$names)


setwd("/home/panos/tmp/mapstat/")   ############ MAGs ############ 
isos=read.table("all.genomic_NOT_imputed.fragment.by.length.ratios_130422")
isos<-isos[!(row.names(isos) %in% "193782"), ]

theix<-subset(sample_info, study=="ours")
season<-as.factor(theix$Sampling.date)
listnames<-rownames(theix)

isos2 <- isos %>% dplyr::select(all_of(listnames)) 
isos2<-t(isos2) # you used not imputed before do the same here
otu.table.filter <- as.data.frame(isos2[,colSums(isos2) > 0])


indval=multipatt(otu.table.filter,season,control = how(nperm=999))
summaryi=indval$sign
summaryi$bh<-p.adjust(summaryi$p.value, method='BH', n = nrow(summaryi))
summaryi$fdr<-p.adjust(summaryi$p.value, method='fdr', n = nrow(summaryi))
summaryi$bonferroni<-p.adjust(summaryi$p.value, method='bonferroni', n = nrow(summaryi))

summary3i <- subset(summaryi, fdr < 0.05)

summary3i$organism=rownames(summary3i)
summary3i$organism<-gsub("(^.*)_k1[0-9]7.*", "\\1", summary3i$organism)
summary3i$organism<-gsub("(^.*).txt.*", "\\1", summary3i$organism)
summary3i$organism<-gsub("(^.*_.*_.*_.*)N[0-9].*", "\\1", summary3i$organism)
summary3i$organism<-gsub("(^.*_.*_.*_.*)W[0-9].*", "\\1", summary3i$organism)
extract<-unique(summary3i$organism)

setwd("/home/panos/tmp")
data<-read.table("phylogeny_metainfo4v2",header=TRUE, row.names = 1, sep=",")
df<-subset(data, data$name %in% extract)
df<-droplevels(df)
df<-df[1:2]
colnames(summary3i)[9]<-c("name")
summary4i<-join(summary3i,df,by="name")
rownames(summary4i)<-rownames(summary3i)

summary4i$contig=rownames(summary4i)
summary4i$contig<-gsub(".txt", ".txt_", summary4i$contig)
summary4i$contig<-gsub("flag.*", "", summary4i$contig)
extract2<-unique(summary4i$contig)


eggNOG.selected<-subset(eggNOG, names %in% extract2)
test<-eggNOG.selected
test$KEGG_ko <- gsub('-',NA,test$KEGG_ko,fixed=TRUE)
test2<-test[complete.cases(test$KEGG_ko),]
colnames(test2)[24]<-c("contig")

Theix<-merge(summary4i,test2,by="contig")

################################################" 

salaheen<-subset(individuality2, study=="salaheen2020")
season<-as.factor(salaheen$diet)
listnames<-rownames(salaheen)

isos2 <- isos %>% dplyr::select(all_of(listnames)) 
isos2<-t(isos2) # you used not imputed before do the same here
otu.table.filter <- as.data.frame(isos2[,colSums(isos2) > 0])


indval=multipatt(otu.table.filter,season,control = how(nperm=999))
summaryi=indval$sign
summaryi$bh<-p.adjust(summaryi$p.value, method='BH', n = nrow(summaryi))
summaryi$fdr<-p.adjust(summaryi$p.value, method='fdr', n = nrow(summaryi))
summaryi$bonferroni<-p.adjust(summaryi$p.value, method='bonferroni', n = nrow(summaryi))

summary3i <- subset(summaryi, fdr < 0.05)

summary3i$organism=rownames(summary3i)
summary3i$organism<-gsub("(^.*)_k1[0-9]7.*", "\\1", summary3i$organism)
summary3i$organism<-gsub("(^.*).txt.*", "\\1", summary3i$organism)
summary3i$organism<-gsub("(^.*_.*_.*_.*)N[0-9].*", "\\1", summary3i$organism)
summary3i$organism<-gsub("(^.*_.*_.*_.*)W[0-9].*", "\\1", summary3i$organism)
extract<-unique(summary3i$organism)

setwd("/home/panos/tmp")
data<-read.table("phylogeny_metainfo4v2",header=TRUE, row.names = 1, sep=",")
df<-subset(data, data$name %in% extract)
df<-droplevels(df)
df<-df[1:2]
colnames(summary3i)[9]<-c("name")
summary4i<-join(summary3i,df,by="name")
rownames(summary4i)<-rownames(summary3i)

summary4i$contig=rownames(summary4i)
summary4i$contig<-gsub(".txt", ".txt_", summary4i$contig)
summary4i$contig<-gsub("flag.*", "", summary4i$contig)
extract2<-unique(summary4i$contig)


eggNOG.selected<-subset(eggNOG, names %in% extract2)
test<-eggNOG.selected
test$KEGG_ko <- gsub('-',NA,test$KEGG_ko,fixed=TRUE)
test2<-test[complete.cases(test$KEGG_ko),]
colnames(test2)[24]<-c("contig")

Salaheen<-merge(summary4i,test2,by="contig")
#################################################



Theix_sep<-subset(subset(Theix, s.sep2020 > 0), stat > 0.5)
Theix_dec<-subset(subset(Theix, s.sep2020 < 1), stat > 0.5)
Salah_milk<-subset(subset(Salaheen, `s.6-milk-solid` < 1), stat > 0.5)
Salah_solid<-subset(subset(Salaheen, `s.6-milk-solid` > 0), stat > 0.5)

Theix %>% tally(s.sep2020 > 0)
Salaheen %>% tally(`s.6-milk-solid` > 0)
(Salaheen %>% tally(`s.6-milk-solid` > 0)*100)/nrow(Salaheen)   # examine percentage of features in each group - this is percentage based on gene counts
(Theix %>% tally(s.sep2020 > 0)*100)/nrow(Theix)





eggNOG$mag<-eggNOG$names
eggNOG$mag<-gsub(".txt_.*","", eggNOG$mag)
eggNOG$mag<-gsub("(^.*)_W[0-9].*", "\\1", eggNOG$mag)
eggNOG$mag<-gsub("(^.*)_N[0-9].*", "\\1", eggNOG$mag)
data<-read.table("phylogeny_metainfo4v2",header=TRUE, row.names = 1, sep = ",")
data<-data[1:3]
colnames(data)[1]<-c("mag")
test7<-merge(eggNOG,data,by="mag",all.x = TRUE)
test8<-subset(test7, !(KEGG_ko == "-"))









#for Theix
mags<-as.data.frame(unique(Theix_sep$original_sample))
mags$cats<-c("sep")
mags2<-as.data.frame(unique(Theix_dec$original_sample))
mags2$cats<-c("dec")
colnames(mags)[1]<-c("original_sample")
colnames(mags2)[1]<-c("original_sample")
maags<-rbind(mags,mags2)



# Theix #
lst<-unique(Theix$original_sample)
test9<-subset(test8, original_sample %in% lst) #that has all KEGGs of the focal MAGs

nrs<-matrix(ncol=2, nrow=length(lst))
colnames(nrs)<-c("labels","completeness")
for(i in 1:length(lst)) {  
  nrs[i,2]<-nrow(subset(Theix_dec, original_sample == lst[i]))/nrow(subset(test9, original_sample == lst[i]))       # this is what you need to check completenesss - gene counts
  nrs[i,1]<-lst[i]} 

Theix_dec_bsetof<-as.data.frame(nrs)
Theix_dec_bsetof$completeness<-as.numeric(Theix_dec_bsetof$completeness)
Theix_dec_bsetof_cntgs<-subset(Theix_dec, original_sample %in% ((subset(Theix_dec_bsetof, completeness > 0.25))$labels))

nrs<-matrix(ncol=2, nrow=length(lst))
colnames(nrs)<-c("labels","completeness")
for(i in 1:length(lst)) {  
  nrs[i,2]<-nrow(subset(Theix_sep, original_sample == lst[i]))/nrow(subset(test9, original_sample == lst[i]))       # this is what you need to check completenesss
  nrs[i,1]<-lst[i]} 

Theix_sep_bsetof<-as.data.frame(nrs)
Theix_sep_bsetof$completeness<-as.numeric(Theix_sep_bsetof$completeness)
Theix_sep_bsetof_cntgs<-subset(Theix_sep, original_sample %in% ((subset(Theix_sep_bsetof, completeness > 0.25))$labels))   #subset MAGs that at least 25% of their genes were diff represented

unique(Theix_sep_bsetof_cntgs$original_sample)
unique(Theix_dec_bsetof_cntgs$original_sample)

median((subset(Theix_sep_bsetof, completeness > 0.25))$completeness)
range((subset(Theix_sep_bsetof, completeness > 0.25))$completeness)

length(unique((subset(Theix_sep, original_sample %in% unique(Theix_sep_bsetof_cntgs$original_sample)))$contig))/length(unique(Theix_sep$contig)) # percentage of contigs






#for Salaheen
mags<-as.data.frame(unique(Salah_milk$original_sample))
mags$cats<-c("milk")
mags2<-as.data.frame(unique(Salah_solid$original_sample))
mags2$cats<-c("solid")
colnames(mags)[1]<-c("original_sample")
colnames(mags2)[1]<-c("original_sample")
maags<-rbind(mags,mags2)


# Salaheen #
lst<-unique(Salaheen$original_sample)
test9<-subset(test8, original_sample %in% lst)

nrs<-matrix(ncol=2, nrow=length(lst))
colnames(nrs)<-c("labels","completeness")
for(i in 1:length(lst)) {  
  nrs[i,2]<-nrow(subset(Salah_milk, original_sample == lst[i]))/nrow(subset(test9, original_sample == lst[i]))       # this is what you need to check completenesss - gene counts
  nrs[i,1]<-lst[i]} 

Salah_milk_bsetof<-as.data.frame(nrs)
Salah_milk_bsetof$completeness<-as.numeric(Salah_milk_bsetof$completeness)
Salah_milk_bsetof_cntgs<-subset(Salah_milk, original_sample %in% ((subset(Salah_milk_bsetof, completeness > 0.1))$labels))

nrs<-matrix(ncol=2, nrow=length(lst))
colnames(nrs)<-c("labels","completeness")
for(i in 1:length(lst)) {  
  nrs[i,2]<-nrow(subset(Salah_solid, original_sample == lst[i]))/nrow(subset(test9, original_sample == lst[i]))       # this is what you need to check completenesss
  nrs[i,1]<-lst[i]} 

Salah_solid_bsetof<-as.data.frame(nrs)
Salah_solid_bsetof$completeness<-as.numeric(Salah_solid_bsetof$completeness)
Salah_solid_bsetof_cntgs<-subset(Salah_solid, original_sample %in% ((subset(Salah_solid_bsetof, completeness > 0.25))$labels))   #subset MAGs that at least 25% of their genes were diff represented

unique(Salah_milk_bsetof_cntgs$original_sample)
unique(Salah_solid_bsetof_cntgs$original_sample)

median((subset(Salah_solid_bsetof, completeness > 0.25))$completeness)
range((subset(Salah_solid_bsetof, completeness > 0.25))$completeness)

length(unique((subset(Salah_solid, original_sample %in% unique(Salah_solid_bsetof_cntgs$original_sample)))$contig))/length(unique(Salah_solid$contig)) # percentage of contigs
