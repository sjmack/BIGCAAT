###Creating a multi-loci data frame with HLA-A and HLA-B

#########################Working with HLA data 
#sets working directory 
setwd("~/Desktop/thesisproject")

#merges combined HLA-A and HLA-B csv files from Bx12 Generator 
hlaABdf<-merge(read.csv("hlaAdownload.csv", header=FALSE, stringsAsFactors=FALSE, skip=3, col.names= c("allelename", "gfes")),
               read.csv("hlaBdownload.csv", header=FALSE, stringsAsFactors=FALSE, skip=3, col.names=c("allelename", "gfes")),
               all=TRUE)

#removes last two rows of combined results 
hlaABdf<-hlaABdf[-seq(nrow(hlaABdf),nrow(hlaABdf)-1),]

#splits allele fields from locus in allele names
#splits allele fields into respective fields 
hlaABfields<-
  strsplit(t(as.data.frame(matrix(unlist(strsplit(hlaABdf[,1],"\\*")), nrow=nrow(hlaABdf), byrow=T), stringsAsFactors=FALSE)) [2,],":")

#finds loci from allele name
#pastes together first two fields to form truncated alleles
#pastes both of these together and into hlaABdf
hlaABdf$match<-paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlaABdf$allelename), paste(sapply(hlaABfields, "[", 1), sapply(hlaABfields, "[", 2), sep=":"), sep="*")


#########################Working with MS relevant data 

#reads in identified MS associated allele values into a data frame for A and B and merges them 
MSABalleles<-merge(read.csv("msAalleles.csv", header=FALSE, stringsAsFactors=FALSE, col.names=c("locus", "MSallele")),
                   read.csv("msBalleles.csv", header=FALSE, stringsAsFactors=FALSE, col.names=c("locus", "MSallele")),
                   all=TRUE)



#matches one gfe to each respective truncated MS allele for A and B 
MSABalleles$GFE<-hlaABdf$gfes[match(paste(MSABalleles$locus, MSABalleles$MSallele, sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlaABdf$allelename), paste(sapply(hlaABfields, "[", 1), sapply(hlaABfields, "[", 2), sep=":"), sep="*")
)]




## Generic Function dealing with multiple files

 MyFunction <- function(fileDir,boolParam=TRUE,textParam="sometext",optionalParam){
  
  combiFile <- multiFileMerge(fileDir)
  View(combiFile[1:10,])
  View(combiFile[5000:5010,])
  
 }
 
 ## From: https://www.r-bloggers.com/merging-multiple-data-files-into-one-data-frame/
 ## Will merge all of the files in mypath/ assuming that they are the same filetype and have common headers 
multiFileMerge <- function(mypath){
   filenames=list.files(path=mypath, full.names=TRUE)
   datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=3,col.names = c("allelename", "gfes"))})
   Reduce(function(x,y) {merge(x,y,all=TRUE)}, datalist) 
   }


