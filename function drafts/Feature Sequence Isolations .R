####Feature Sequence Isolation for GFEs 
##By: Livia Tran 
#V 1.3

###loads pre-made RDA files 
###framework.rda is a list consisting of feature identification for each field of a GFE
###atlas.rda is a dataframe that gives the position for selected features in a GFE
load("/Users/liviatran/ltmasterscoding/framework.rda")
load("/Users/liviatran/ltmasterscoding/atlas.rda")

#loads necessary library for usage of str_replace_all
library(stringr)

##function for reading multiple HLA files in, but keeping them as separate dataframes
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/HLA", c("allelename", "gfe"), skip=3, clip=1)

#splits GFE notations -- first at "w", then at "-"
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

#loops through rows in atlas -- takes each row in "locus", and splits the coinciding allele,
#forming a dataframe
#inputs into a blank splitDFs list
#field names, previously defined by framework, corresponds to HLA locus -- inserts 
#appropriate names to splitDFs
splitDFs<-list()
for(i in 1:nrow(atlas)){
  splitDFs[[atlas$locus[i]]]<-data.frame(matrix(unlist(hlasplit[i]), byrow=T, nrow=length(hlasplit[[i]])), stringsAsFactors = FALSE)
  colnames(splitDFs[[atlas$locus[i]]])<-framework[[i]]}


#function to make zeroed out templates based on each loci - then fills in with appropriate information
#based on coordinates of the atlas 
featureselect<-function(x, y, GFEsplitDF){
  dataframe <- GFEsplitDF
  dataframe[,1:ncol(dataframe)] <- as.data.frame(rbind(rep(0,length(dataframe))))
  dataframe[,unlist(atlas[x,y])] <- GFEsplitDF[unlist(atlas[x,y])]
  return(dataframe)}

#list.names pastes column names of sequence features, followed by "list"
list.names<-paste(colnames(atlas[,2:9]), "list", sep="")

#makes an empty list of 8 total features, names of elements previously defined by
#list.names
seqfeatureslist <- sapply(list.names,function(x) NULL)

#inputs UTR isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[1]][[atlas$locus[i]]]<-featureselect(i ,2, splitDFs[[i]])}
#inputs core exon isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[2]][[atlas$locus[i]]]<-featureselect(i ,3, splitDFs[[i]])}
#inputs intercytoplasmic exon isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[3]][[atlas$locus[i]]]<-featureselect(i ,4, splitDFs[[i]])}
#inputs leader peptide exon isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[4]][[atlas$locus[i]]]<-featureselect(i ,5, splitDFs[[i]])}
#inputs c-domain exon isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[5]][[atlas$locus[i]]]<-featureselect(i ,6, splitDFs[[i]])}
#inputs intron isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[6]][[atlas$locus[i]]]<-featureselect(i ,7, splitDFs[[i]])}
#inputs transmembrane exon isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[7]][[atlas$locus[i]]]<-featureselect(i ,8, splitDFs[[i]])}
#inputs non-core exon isolations for all loci in list format 
for (i in 1:nrow(atlas)){
  seqfeatureslist[[8]][[atlas$locus[i]]]<-featureselect(i ,9, splitDFs[[i]])}

#function for pasting all GFE fields together
#adds an additional row to all loci for each sequence feature
#inputs a *00:00 allele, with a maxed out GFE for each field 
GFEpaster<-function(seqfeaturelist, seqfeatureposition, lociDF, BSGdf, BSGdfposition, zerolocus){
  seqfeaturelist[[seqfeatureposition]][nrow(seqfeaturelist[[seqfeatureposition]]) + 1,] = max(as.numeric(as.character(unlist(lociDF[,2:ncol(lociDF)]))))+1
  cols<- 1:length(seqfeaturelist[[seqfeatureposition]])
  GFEs<-str_replace_all(paste(gsub("[*].*", "", unlist(BSGdf[[BSGdfposition]][1])),paste("w", apply(seqfeaturelist[[seqfeatureposition]][,cols] , 1 , paste, collapse="-"))), " ", "")
  BSGdf[[BSGdfposition]][nrow(BSGdf[[BSGdfposition]]) + 1,] = list(zerolocus, NA)
  featureselectDF<-data.frame(cbind(BSGdf[[BSGdfposition]][[1]], GFEs))
  colnames(featureselectDF)[1] <- colnames(BSGdf[[BSGdfposition]][,])[1]
  return(featureselectDF)}

##isolated.names pastes column names of sequence features, followed by "isolated"
isolated.names<-paste(colnames(atlas[,2:9]), "isolated", sep="")

#makes an empty list of 8 total features, names of elements previously defined by
#isolated.names
isolatedlist <- sapply(isolated.names,function(x) NULL)

#for loop for pasting all sequence feature fields together
#inputs into isolatedlist
for(i in 1:nrow(atlas)){
  for(k in 1:length(seqfeatureslist)){
  isolatedlist[[k]][[atlas$locus[i]]]<-GFEpaster(seqfeatureslist[[k]], i, splitDFs[[i]], hladf, i, str_replace_all(paste(atlas$locus[i], "*00:00"), " ", ""))}}





