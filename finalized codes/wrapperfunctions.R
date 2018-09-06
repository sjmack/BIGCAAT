###functions 

require(stringr)
require(BIGDAWG)

#for file merging
filemerge <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<-Reduce(function(x,y) {merge(x,y,all=TRUE)}, datalist)
  datalist<-head(datalist,nrow(datalist)-(clip*length(filenames)))
  return(datalist)}

#for cwd pairing 
CWDrestriction <- function(allelelistFile,cwdFile,gfeDirPath) {
  hlaacc<-read.csv(allelelistFile, header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")
  hladf<-filemerge(gfeDirPath, c("allelename", "gfe"), skip=3, clip=1)
  hladf$alleleID<-hlaacc$AlleleID[match(hladf$allelename, paste("HLA",hlaacc$Allele,sep="-"))] 
  hladf$CWD<-ifelse(hladf$alleleID %in% (read.table(cwdFile,sep="\t",header = TRUE,stringsAsFactors = FALSE,skip = 1)$IMGT.HLA.Accession.Number), "CWD", "NON-CWD")
  return(hladf)}

#for MS allele GFE matching -- hlamerged previously defined by CWDrestriction
gfematch <- function(allelefiles,cwddata,gfepath, experimentaldata){
  hlamerged<-CWDrestriction(allelelistFile=allelefiles, cwdFile=cwddata, gfeDirPath=gfepath)
  msmerged<-filemerge(experimentaldata, c("locus", "MSallele"), 0, 0)
  hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")
  msmerged$gfeCWD<-ifelse(hlamerged$CWD[match(str_replace_all(string=(paste(msmerged$locus, msmerged$MSallele, sep="*")), pattern=" ", repl=""), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="CWD", 
                          hlamerged$gfe[match(str_replace_all(string=(paste(msmerged$locus, msmerged$MSallele, sep="*")), pattern=" ", repl=""), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                          )], ifelse(hlamerged$CWD[match(str_replace_all(string=(paste(msmerged$locus, msmerged$MSallele, sep="*")), pattern=" ", repl=""), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="NON-CWD", 
                                     hlamerged$gfe[match(str_replace_all(string=(paste(msmerged$locus, msmerged$MSallele, sep="*")), pattern=" ", repl=""), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                                     )], "")) 
  return(msmerged)}

#for multi-loci GFE genotype matching 
GenotypeGFEmatch <- function(allelefiles,cwddata,gfepath, genotypedatafiles){
  hlamerged<-CWDrestriction(allelefiles,cwddata,gfepath)
  hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")
  genodata<-filemerge(genotypedatafiles, "genodata" , 0, 0)
  genofields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(genodata[,1],"\\*")), nrow=nrow(genodata), byrow=T), stringsAsFactors=FALSE)) [2,],":")
  genodata$gfeCWD<-ifelse(hlamerged$CWD[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="CWD", 
                          hlamerged$gfe[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                          )], ifelse(hlamerged$CWD[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="NON-CWD", 
                                     hlamerged$gfe[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                                     )], "")) 
  return(genodata)}

#BDgenotypeconversion
BDgenotypeconversion<-function(gfeDirpath, genotypedata, allelefiles, cwdfiles){
if (is.data.frame(genotypedata)) {bigdawghladata <- genotypedata} else 
{switch(tolower(substr(genotypedata,nchar(genotypedata)-2,nchar(genotypedata))),
                  "csv" = bigdawghladata <- read.csv(genotypedata, stringsAsFactors = F, header=T, check.names = F), 
                  "txt" = bigdawghladata <- read.table(genotypedata, stringsAsFactors = F, header=T, sep="\t", check.names = FALSE), 
                  "tsv" = bigdawghladata <- read.table(genotypedata, stringsAsFactors = F, header=T,sep="\t",check.names = FALSE), 
                  print("Error: Unrecognized filename suffix."))}
hlamerged<-CWDrestriction(allelefiles, cwdfiles, gfeDirpath)
hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")
hlalist<-list()            
  for (i in 3:ncol(bigdawghladata)){
    hlalist[[i]]<-strsplit(bigdawghladata[,i],":")} 
   for(i in 3:ncol(bigdawghladata)){
    bigdawghladata[i]<-paste(sapply(hlalist[[i]], "[", 1), sapply(hlalist[[i]], "[", 2), sep=":")
        bigdawghladata[,i]<-paste(colnames(bigdawghladata[i]),bigdawghladata[,i],sep="*")}
   for (i in 3:ncol(bigdawghladata)){
      bigdawghladata[i]<-hlamerged$gfe[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]}
return(bigdawghladata)
  }


setwd("/Users/liviatran/Desktop/ltmasterscoding")
View(BDgenotypeconversion("/Users/liviatran/Desktop/ltmasterscoding/HLA", "/Users/liviatran/Desktop/ltmasterscoding/hladata.txt", "Allelelist.3310.txt","cwd200_alleles.txt"))
View(gfematch("Allelelist.3310.txt","cwd200_alleles.txt","/Users/liviatran/Desktop/ltmasterscoding/HLA","/Users/liviatran/Desktop/ltmasterscoding/MS"))
