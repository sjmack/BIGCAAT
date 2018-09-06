#sets working directory 
setwd("/Users/liviatran/Desktop/ltmasterscoding")

#needed packages 
require(BIGDAWG)

#loads pre-existing HLA data from BIGDAWG
hladata<-HLA_data

#merges hla files
hlamerged<-CWDrestriction("Allelelist.3310.txt","cwd200_alleles.txt","/Users/liviatran/Desktop/ltmasterscoding/HLA")

#splits into 4 elements
hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")

#for loop for splitting each column -- makes into a large list 
hlalist<-list()            
for (i in 3:ncol(hladata)){
  hlalist[[i]]<-strsplit(hladata[,i],":")}            


#for loop for pasting first two fields together to make truncated allele
for(i in 3:ncol(hladata)){
hladata[i]<-paste(sapply(hlalist[[i]], "[", 1), sapply(hlalist[[i]], "[", 2), sep=":")
hladata[,i]<-paste(colnames(hladata[i]),hladata[,i],sep="*")}

#for loop for determining which genotype alleles are CWD and whcih are not
CWD<-cbind(hladata)
for(i in 3:ncol(hladata)){
CWD[i] <- hlamerged$CWD[match(hladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]}

#for loop for matching GFEs to CWD alleles - if not a CWD allele, dubs "non-CWD"
for(i in 3:ncol(hladata)){
hladata[i]<-ifelse(CWD[,i]=="CWD", hlamerged$gfe[match(hladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))], ifelse(CWD[,i]=="NON-CWD", hlamerged$gfe[match(hladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))], "NA")
)}
           

#filemerge function 
filemerge <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<-Reduce(function(x,y) {merge(x,y,all=TRUE)}, datalist)
  datalist<-head(datalist,nrow(datalist)-(clip*length(filenames)))
  return(datalist)}





