setwd("~/Desktop/ltmasterscoding")

require(stringr)

#reads in appropriate allele text file
#separates the first row into Allele and AlleleID -- makes into a dataframe
#renames the columns based on the first row (i.e. what was originally in the text file)
#deletes first row that was made into a header 
hlaacc<-read.csv("Allelelist.3310.txt", header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")

#reads in cwd tab-delimited data
#renames headers from those found in the first row
#deletes first row that was made into a header 
cwddata <- read.table("cwd200_alleles.txt",sep="\t",header = TRUE,stringsAsFactors = FALSE,skip = 1)

#merges hla files based on filemerge function 
hladf<-filemerge("/Users/liviatran/Desktop/ltmasterscoding/HLA", c("allelename", "gfe"), skip=3, clip=1)

#adds "HLA" in front of alleles in hlaacc
#matches allelenames to those found in hlaacc to find alleleIDs for each allele name
#further subsets hladf to those found in the CWD database 
hladf<-subset(hladf, hlaacc$AlleleID[match(hladf$allelename, paste("HLA",hlaacc$Allele,sep="-"))] %in% cwddata$IMGT.HLA.Accession.Number)



#filemerge function 
filemerge <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<-Reduce(function(x,y) {merge(x,y,all=TRUE)}, datalist)
  datalist<-head(datalist,nrow(datalist)-(clip*length(filenames)))
  return(datalist)}
                        