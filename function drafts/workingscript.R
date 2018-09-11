########UTR SEQUENCE FEATURE SELECT

##required packages
require(BIGDAWG)
require(stringr)

###necessary functions, with descriptions, for this script are found at the bottom of the script###

##see atlas.r, which contains positions for hand-picked sequence features
load("atlas.rda")

#BDgenotypeconversion was run on HLA_data in the BIGDAWG package for this script -- 
#samplebdHLA.csv was a downloaded CSV version of HLA_data
hla_data<-BDgenotypeconversion("/Users/liviatran/Desktop/ltmasterscoding/HLA", "samplebdHLA.csv", "Allelelist.3310.txt", CWDverify())

##creates a new list of hla_data without the NAs (i.e. missing values)
hla_omitted<-list()
for(i in 3:ncol(hla_data)){
  hla_omitted[[i]]<-na.omit(hla_data[,i])}

##hlaA
hlaA<-combine.lists(hla_omitted[[3]], hla_omitted[[4]])
hlaAdf<-data.frame(Reduce(rbind, strsplit(t(as.data.frame(matrix(unlist(strsplit(hlaA,"w")), nrow=length(hlaA), byrow=T), stringsAsFactors=FALSE)) [2,], "-")
))
colnames(hlaAdf)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                  "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR") 


#hlaDRB1
hlaDRB1<-combine.lists(hla_omitted[[5]], hla_omitted[[6]])
hlaDRB1df<-data.frame(Reduce(rbind, strsplit(t(as.data.frame(matrix(unlist(strsplit(hlaDRB1,"w")), nrow=length(hlaDRB1), byrow=T), stringsAsFactors=FALSE)) [2,], "-")
))
colnames(hlaDRB1df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")



##prepares for building a list for UTR selected features for each loci
loci<-c("HLA-A", "HLA-DRB1", "HLA-DQB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5")
UTRlist<-list()


A_UTR_template<-hlaAdf
A_UTR_template[,1:ncol(A_UTR_template)] <- as.data.frame(rbind(rep(0,17)))

DRB1_UTR_template<-hlaDRB1df
DRB1_UTR_template[,1:ncol(DRB1_UTR_template)] <- as.data.frame(rbind(rep(0,13)))

UTRlist[[loci[1]]]<-atlasselect(A_UTR_template, 1,2, hlaAdf) 
UTRlist[[loci[2]]]<-atlasselect(DRB1_UTR_template, 8,2, hlaDRB1df) 


View(UTRlist[2])




############## NEEDED FUNCTIONS 
##function to select specific sequence features based on their position, which are documented by atlas.R
##GFEsplitDF denotes a loci specific dataframe, with each GFE field having its own column 
atlasselect<-function(seqfeature, x, y, GFEsplitDF){
  seqfeature[,unlist(atlas[x,y])] <- GFEsplitDF[unlist(atlas[x,y])]
  return(seqfeature)}

#for splitting each specific element of a list 
listSplit<-function(somelist, element){
  splitelement<-strsplit(t(as.data.frame(matrix(unlist(strsplit(somelist[[element]],"w")), nrow=length(somelist[[element]]), byrow=T), stringsAsFactors=FALSE)) [2,], "-")
  return(splitelement)
}

#for file merging
filemerge <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<-Reduce(function(x,y) {merge(x,y,all=TRUE)}, datalist)
  datalist<-head(datalist,nrow(datalist)-(clip*length(filenames)))
  return(datalist)}

#for common, well-documented allele identification 
cwdID <- function(allelelistFile,gfeDirPath) {
  hlaacc<-read.csv(allelelistFile, header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")
  hladf<-filemerge(gfeDirPath, c("allelename", "gfe"), skip=3, clip=1)
  hladf$alleleID<-hlaacc$AlleleID[match(hladf$allelename, paste("HLA",hlaacc$Allele,sep="-"))] 
  cwdalleles<-CWDverify()
  hladf$CWD<-ifelse(hladf$alleleID %in% cwdalleles$Accession, "CWD", "NON-CWD")
  return(hladf)}

#for converting multi-loci data to its GFE counterpart
BDgenotypeconversion<-function(gfeDirpath, genotypedata, allelefiles, cwdfiles){
  if (is.data.frame(genotypedata)) {bigdawghladata <- genotypedata} else 
  {switch(tolower(substr(genotypedata,nchar(genotypedata)-2,nchar(genotypedata))),
          "csv" = bigdawghladata <- read.csv(genotypedata, stringsAsFactors = F, header=T, check.names = F), 
          "txt" = bigdawghladata <- read.table(genotypedata, stringsAsFactors = F, header=T, sep="\t", check.names = FALSE), 
          "tsv" = bigdawghladata <- read.table(genotypedata, stringsAsFactors = F, header=T,sep="\t",check.names = FALSE), 
          print("Error: Unrecognized filename suffix."))}
  hlamerged<-cwdID(allelefiles, gfeDirpath)
  hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")
  hlalist<-list()             
  for (i in 3:ncol(bigdawghladata)){
    hlalist[[i]]<-strsplit(bigdawghladata[,i],":")} 
  for(i in 3:ncol(bigdawghladata)){
    bigdawghladata[i]<-paste(sapply(hlalist[[i]], "[", 1), sapply(hlalist[[i]], "[", 2), sep=":")
    bigdawghladata[,i]<-paste(colnames(bigdawghladata[i]),bigdawghladata[,i],sep="*")}
  for (i in 3:ncol(bigdawghladata)){
    bigdawghladata[i]<-ifelse(hlamerged$CWD[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="CWD", 
                              hlamerged$gfe[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                              )], ifelse(hlamerged$CWD[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="NON-CWD", 
                                         hlamerged$gfe[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1",  hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                                         )], ""))} 
  
  return(bigdawghladata)}


##CWDverify() function, authored by Steve J. Mack -- this function pulls the most recent CWD files and compares them to
#alleles that my have changed in name or have been deleted -- outputs CWD accession numbers and alleles, updated with
#any possible name changes or deletions 
CWDverify <- function(){
  require(data.table)
  
  ## Pull down the CWD catalogue
  CWD <- list()
  CWD$data <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",skip = 1,stringsAsFactors = FALSE,select = c(2,3))
  CWD$version <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",nrows = 1,stringsAsFactors = FALSE,select=1)
  
  ## Pull down the hla_nom.txt, Deleted_alleles.txt and allelelist.txt files to create a table of v3.0.0+ deleted alleles, their ACCs,their replacements, and their ACCs
  deletedHLA <- list()
  # Temporarily store the entire hla_nom.txt in $version
  deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",skip=6, stringsAsFactors = FALSE,sep = ";", col.names = c("Locus","AlleleName","NewName","Event"),select = c(1,2,5,6))
  ## Exclude entries without allele name changes
  deletedHLA$data <- deletedHLA$version[deletedHLA$version$NewName !="",]
  # Exclude pre-db release 3.0.0 alleles
  deletedHLA$data <- deletedHLA$data[grep(":",deletedHLA$data$AlleleName,fixed=TRUE),]
  
  ## Process and extract the accession numbers from the Deleted_alleles.txt file, temporarily stored in $version
  deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,skip = 7,sep=",",header=TRUE,fill=TRUE)
  ## Below to account for one extra comma in line 106 (hopefully, can be deleted in a future release)
  if(ncol(deletedHLA$version)==4) {deletedHLA$version$Description[98] <- paste(deletedHLA$version$Description[98],deletedHLA$version$V4[98],sep=" ")
  deletedHLA$version <- deletedHLA$version[,1:3] }
  # Store the pertinent accession numbers in the data element
  deletedHLA$data$origAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$AlleleName,sep=""),deletedHLA$version$Allele)]
  # Temporarily store the allelelist.txt file in $version 
  deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",skip=6, stringsAsFactors = FALSE,sep = ",", header=TRUE)
  deletedHLA$data$newAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep=""),deletedHLA$version$Allele)]
  # overwrite the Deleted_alelles.txt files with the version information
  deletedHLA$version <- cbind(fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE),fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE), fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",nrows=5, stringsAsFactors = FALSE,sep = "?", header=TRUE))
  
  ## Match accession numbers in CWD to the Accession numbers in the deleted alleles. 
  changeCWD <- match(CWD$data$`IMGT/HLA Accession Number`,deletedHLA$data$origAccession)
  # Create full allele names for the new names
  deletedHLA$data$NewName <- paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep="")
  CWD$data[!is.na(changeCWD),] <- cbind(deletedHLA$data[changeCWD[!is.na(changeCWD)],6],deletedHLA$data[changeCWD[!is.na(changeCWD)],3])
  
  # Rename the columns of the verified CWD table
  colnames(CWD$data) <- c("Accession","AlleleName")
  
  CWD$data
}


