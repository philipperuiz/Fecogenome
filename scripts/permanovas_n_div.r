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

#all files can be downloaded at: https://doi.org/10.15454/VWZBUF
setwd("/home/panos/tmp/")

sample_info<-read.table("categories_samples2", sep="\t", header=TRUE, row.names = 1,quote="", fill=FALSE)
individuality=subset(sample_info, sample.type=="individual")
#individuality<-individuality[complete.cases(individuality$diet), ]

theix<-subset(sample_info, study=="ours")
individuality2<-individuality[complete.cases(individuality$diet),]
individuality=droplevels(individuality)
individuality2=droplevels(individuality2)
salaheen<-subset(individuality2, study=="salaheen2020")


individuality2$location<-as.factor(individuality2$location)
individuality2$Animal<-as.factor(individuality2$Animal)
individuality2$Sampling.date<-as.factor(individuality2$Sampling.date)
individuality2$diet<-as.factor(individuality2$diet)
individuality2$study<-as.factor(individuality2$study)
length(levels(individuality2$study))
individuality2$merged<-paste(individuality2$study,individuality2$Animal)
individuality2$merged<-as.factor(individuality2$merged)
levels(individuality2$merged)

individuality2$Animal2<-individuality2$Animal
individuality2$Animal2<-gsub("lees than 2yo steer", "M", individuality2$Animal2)
individuality2$Animal2<-as.factor(individuality2$Animal2)
levels(individuality2$Animal2)






setwd("/home/panos/tmp/mapstat")   ############ MAGs ############ 
isos=read.table("all.genomic.corrected_n_imputed.fragment.by.length.ratios_130422")


isos2<-isos[!(row.names(isos) %in% "193782"), ]
colnames(isos2)[109]<-c("N15")
listnames=rownames(individuality2)
isos3 <- isos2 %>% dplyr::select(all_of(listnames))

isos5=read.table("all.genomic_NOT_imputed.fragment.by.length.ratios_130422")    # CORRECTD BUT NOT IMPUTED
colnames(isos5)[109]<-c("N15")
isos6<-isos5[!(row.names(isos5) %in% "193782"), ]
isos7 <- isos6 %>% dplyr::select(all_of(listnames))
isos7=t(isos7)
isos7=as.data.frame(isos7)
filter <- isos7[,colSums(isos7) > 0]
names<-colnames(filter)
otu.table.filter <- subset(isos3, rownames(isos3) %in% names)
otu.table.filter<-droplevels(otu.table.filter)


otu.table.filter.clr<-as.data.frame(clr(t(otu.table.filter)))


############ PERMANOVASs ############" 

set.seed(123)
adonis(data=individuality2, tesssst ~ diet, permutations = 999, method = "euclidean" )
tbl4<-adonis(data=individuality2, otu.table.filter.clr ~ study:Animal2, permutations = 999, method = "euclidean" ) # should be the same
tbl4$aov.tab


#theix
ours<-rownames(subset(individuality, study=="ours"))
otu.table.filter<-as.data.frame(t(otu.table.filter))
otu.table.filter.filter<-subset(otu.table.filter, rownames(otu.table.filter) %in% ours)
otu.table.filter.filter<-droplevels(otu.table.filter.filter)
otu.table.filter.filter.clr<-as.data.frame(clr(otu.table.filter.filter))

data<-subset(individuality2, study=="ours")
data<-droplevels(data)

data$host<-as.factor(data$host)
set.seed(9)
tbl2<-adonis(data=data, otu.table.filter.filter.clr ~ Sampling.date+host, permutations = 999, method = "euclidean" )
tbl2$aov.tab


#salaheen
data<-subset(individuality2, study=="salaheen2020")
data<-droplevels(data)
otu.table.filter.filter<-subset(otu.table.filter, rownames(otu.table.filter) %in% rownames(data))
otu.table.filter.filter<-droplevels(otu.table.filter.filter)
otu.table.filter.filter.clr<-as.data.frame(clr(otu.table.filter.filter))
set.seed(0)
tbl<-adonis(data=data, otu.table.filter.filter.clr ~ diet+Farm, permutations = 999, method = "euclidean" )
tbl$aov.tab
salah<-tbl$aov.tab





######################## alpha diversity for MAGs ######################## 
otus<-subset(isos6, rownames(isos6) %in% names)   # using the corrected and not imputed
otus<-t(otus)
otus<-subset(otus, rownames(otus) %in% listnames)
otus<-as.data.frame(otus)
otus<-droplevels(otus)

shannon<-as.data.frame(diversity(otus), index = "shannon")
Sobs<-as.data.frame(log(colSums(t(otus) != 0)))

isossss<-cbind(shannon,Sobs)
isossss$sample<-rownames(isossss)
colnames(isossss)<-c("shannon","logSobs","sample")

individuality2$sample<-rownames(individuality2)
isossss3<-merge(isossss,individuality2,by="sample")

#theix
temp<-subset(isossss3, study=="ours")
temp<-droplevels(temp)

kruskal.test(logSobs ~ Sampling.date, data = temp)
kruskal.test(shannon ~ Sampling.date, data = temp)

#salaheen
temp<-subset(isossss3, study=="salaheen2020")
temp<-droplevels(temp)

kruskal.test(logSobs ~ diet, data = temp)
kruskal.test(shannon ~ diet, data = temp)








####################### alpha diversity for virulence features ####################### 

setwd("/home/panos/tmp/mapstat")

new3<-read.table("VFdb.VFcorrected.NOT_imputed.fragment.by.length.ratios_130422")
colnames(new3)<-gsub("_not_host_VFdb.mapstat4","", colnames(new3))
new3<-new3[setdiff(rownames(new3),"colSums(new3)"),]      # remove useless last row # CHECK

new3$organism=rownames(new3)
new3$pathogene=rownames(new3)

new3$organism<-gsub("\\]", "£", gsub("\\[", "@", new3$organism))
new3$pathogene<-gsub("\\]", "£", gsub("\\[", "@", new3$pathogene))
new3$pathogene<-gsub("(^.*)(@.*£)(@.*£)", "\\2", new3$pathogene)
new3$organism<-gsub("(^.*)(@.*£)(@.*£)", "\\3", new3$organism)

new3$organism<-gsub("@", "", new3$organism)
new3$organism<-gsub("£", "", new3$organism)
new3$pathogene<-gsub("@", "", new3$pathogene)
new3$pathogene<-gsub("£", "", new3$pathogene)


new3$querry<-rownames(new3)





list1<-c("flaA", "cadF", "iam", "cdt", "virB11", "wlaN")
list2<-c("Campylobacter")
newToxin1<-new3[grep(paste0(list1, collapse = "|"), new3$querry), ]
newToxin2<-new3[grep(paste0(list2, collapse = "|"), new3$querry), ] 
newToxin1.2<-newToxin2[grep(paste0(list1, collapse = "|"), newToxin2$querry), ] 
campy<-newToxin1.2

list1<-c("eae", "bfpA", "pic", "papG", "afaB", "afaC", "aggR", "ipaH", "daaD", "aaiC", "aatA", "stable enterotoxin")
list2<-c("Escherichia","Shigella")
newToxin1<-new3[grep(paste0(list1, collapse = "|"), new3$querry), ]
newToxin2<-new3[grep(paste0(list2, collapse = "|"), new3$querry), ] 
newToxin1.2<-newToxin2[grep(paste0(list1, collapse = "|"), newToxin2$querry), ] 
Esche_shig<-newToxin1.2

list1<-c("cpA", "cpB", "cpE","cpa", "cpb", "cpe")
list2<-c("Clostridium")
newToxin1<-new3[grep(paste0(list1, collapse = "|"), new3$querry), ]
newToxin2<-new3[grep(paste0(list2, collapse = "|"), new3$querry), ] 
newToxin1.2<-newToxin2[grep(paste0(list1, collapse = "|"), newToxin2$querry), ] 
clostr<-newToxin1.2

list1<-c("invA", "fimA", "stn", "spvC", "spvR")
list2<-c("Salmonella")
newToxin1<-new3[grep(paste0(list1, collapse = "|"), new3$querry), ]
newToxin2<-new3[grep(paste0(list2, collapse = "|"), new3$querry), ] 
newToxin1.2<-newToxin2[grep(paste0(list1, collapse = "|"), newToxin2$querry), ] 
salmo<-newToxin1.2

list1<-c("internalin", "hly", "plcB", "plcA", "ActA", "iap")
list2<-c("monocytogenes")
newToxin1<-new3[grep(paste0(list1, collapse = "|"), new3$querry), ]
newToxin2<-new3[grep(paste0(list2, collapse = "|"), new3$querry), ] 
newToxin1.2<-newToxin2[grep(paste0(list1, collapse = "|"), newToxin2$querry), ] 
liste<-newToxin1.2

temp2<-rbind(campy,clostr,Esche_shig,liste,salmo)



shannon<-as.data.frame(diversity(t(temp2[1:436]), index = "shannon"))
Sobs<-as.data.frame(log(colSums(temp2[1:436] != 0)))
colnames(Sobs)<-c("val")
Sobs$val[which(!is.finite(Sobs$val))] <- 0 # replace "Inf" with "0s" glitch happens cause these samples have only 0s so diversity is 0



isossss<-cbind(shannon,Sobs)
isossss$sample<-rownames(isossss)
colnames(isossss)<-c("shannon","logSobs","sample")


isossss2<-subset(isossss, rownames(isossss) %in% rownames(individuality))
isossss3<-cbind(isossss2,individuality)



#theix
temp<-subset(isossss3, study=="ours")
temp<-droplevels(temp)

kruskal.test(logSobs ~ Sampling.date, data = temp)
kruskal.test(shannon ~ Sampling.date, data = temp)

#salaheen
temp<-subset(isossss3, study=="salaheen2020")
temp<-droplevels(temp)

kruskal.test(logSobs ~ diet, data = temp)
kruskal.test(shannon ~ diet, data = temp)









