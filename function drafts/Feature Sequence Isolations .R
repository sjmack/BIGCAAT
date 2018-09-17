####Feature Sequence Isolation for GFEs 
##By: Livia Tran 
#V 1.4
#September 17, 2018


############NECESSARY FUNCTIONS 
#CWDverify, by SJ MACK -- compares CWD catalogue to changed and deleted files
#outputs accurate CWD files 
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

#for identifying common, well-documented alleles in HLA data 
cwdID <- function(allelelistFile,gfefiles) {
  hlaacc<-read.csv(allelelistFile, header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")
  hladf<-gfefiles
  hladf$alleleID<-hlaacc$AlleleID[match(hladf$allelename, paste("HLA",hlaacc$Allele,sep="-"))] 
  hladf$alleleID<-ifelse(is.na(hladf$alleleID)=="TRUE", paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hladf$allelename), "x", sep=""), hladf$alleleID)
  cwdalleles<-CWDverify()
  hladf$CWD<-ifelse(hladf$alleleID %in% cwdalleles$Accession, "CWD", "NON-CWD")
  return(hladf)}

View(cwdID("/Users/liviatran/Desktop/ltmasterscoding/Allelelist.3310.txt", mergedlist[[1]]))
##function for reading multiple HLA files in, but keeping them as separate dataframes
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}

#function to make zeroed out templates based on each locus - then fills in with appropriate information
#based on coordinates of the atlas 
featureselect<-function(x, y, GFEsplitDF){
  dataframe <- GFEsplitDF
  dataframe[,1:ncol(dataframe)] <- as.data.frame(rbind(rep(0,length(dataframe))))
  dataframe[,unlist(atlas[x,y])] <- GFEsplitDF[unlist(atlas[x,y])]
  return(dataframe)}

#function for pasting all GFE fields together
#adds an additional row to all loci for each sequence feature
#inputs a *00:00 allele, with a maxed out GFE for each field 
GFEpaster<-function(seqfeaturelist, seqfeatureposition, lociDF, BSGdf, BSGdfposition, zerolocus){
  seqfeaturelist[[seqfeatureposition]][nrow(seqfeaturelist[[seqfeatureposition]]) + 1,] = max(as.numeric(as.character(unlist(lociDF[,2:ncol(lociDF)]))))+1
  cols<- 1:length(seqfeaturelist[[seqfeatureposition]])
  GFEs<-str_replace_all(paste(gsub("[*].*", "", unlist(BSGdf[[BSGdfposition]][1])),paste("w", apply(seqfeaturelist[[seqfeatureposition]][,cols] , 1 , paste, collapse="-"))), " ", "")
  BSGdf[[BSGdfposition]][nrow(BSGdf[[BSGdfposition]]) + 1,] = list(zerolocus, NA)
  featureselectDF<-cbind.data.frame(BSGdf[[BSGdfposition]][[1]], GFEs, stringsAsFactors=FALSE)
  colnames(featureselectDF)[1] <- colnames(BSGdf[[BSGdfposition]][,])[1]
  return(featureselectDF)}



###########BEGIN SCRIPT 

###loads pre-made RDA files 
###framework.rda is a list consisting of feature identification for each field of a GFE
###atlas.rda is a dataframe that gives the position for selected features in a GFE
load("/Users/liviatran/ltmasterscoding/framework.rda")
load("/Users/liviatran/ltmasterscoding/atlas.rda")

#loads necessary libraries
require(stringr)
require(BIGDAWG)

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


#makes an empty list of 8 total features, names of elements defined by
#columnnames of the atlas
seqfeatureslist <- sapply(colnames(atlas[,2:9]),function(x) NULL)

#for loop for inputting sequence feature isolations for all loci in list format 
for (k in 1:length(seqfeatureslist)) {
  for (i in 1:nrow(atlas)){
    seqfeatureslist[[k]][[atlas$locus[i]]]<-featureselect(i ,k+1, splitDFs[[i]])}
}

#makes an empty list of 8 total features, names of elements defined by atlas column names
isolatedlist <- sapply(colnames(atlas[,2:9]),function(x) NULL)

#for loop for pasting all sequence feature fields together
#inputs into isolatedlist
for(i in 1:nrow(atlas)){
  for(k in 1:length(seqfeatureslist)){
    isolatedlist[[k]][[atlas$locus[i]]]<-GFEpaster(seqfeatureslist[[k]], i, splitDFs[[i]], hladf, i, as.character(str_replace_all(paste(atlas$locus[i], "*00:00"), " ", "")))}}

##makes an empty list of 8 total features, names of elements defined by atlas column names
mergedlist<-sapply(colnames(atlas[,2:9]),function(x) NULL)

#recursively merges all dataframes for each sequence feature together from isolated list
#inputs into merged list 
for(i in 1:length(isolatedlist)){
mergedlist[[i]]<-as.data.frame(Reduce(function(x,y) {merge(x,y,all=TRUE)}, isolatedlist[[i]]))}   


BDgenotypeconversion<-function(genotypedata, allelefiles, gfefiles){
  proceed <- TRUE
  if (is.data.frame(genotypedata)) {bigdawghladata <- genotypedata} else 
  {switch(tolower(substr(genotypedata,nchar(genotypedata)-2,nchar(genotypedata))),
          "csv" = bigdawghladata <- read.csv(genotypedata, stringsAsFactors = F, header=T, check.names = F), 
          "txt" = bigdawghladata <- read.table(genotypedata, stringsAsFactors = F, header=T, sep="\t", check.names = FALSE), 
          "tsv" = bigdawghladata <- read.table(genotypedata, stringsAsFactors = F, header=T,sep="\t",check.names = FALSE), 
          proceed <- FALSE)}
  if(proceed) {
    hlamerged<-cwdID(allelefiles, gfefiles)
    hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")
    hlalist<-list()            
    for (i in 3:ncol(bigdawghladata)){
      hlalist[[i]]<-strsplit(bigdawghladata[,i],":")} 
    for(i in 3:ncol(bigdawghladata)){
      bigdawghladata[i]<-paste(sapply(hlalist[[i]], "[", 1), sapply(hlalist[[i]], "[", 2), sep=":")
      bigdawghladata[,i]<-paste(colnames(bigdawghladata[i]),bigdawghladata[,i],sep="*")}
    for (i in 3:ncol(bigdawghladata)){
      bigdawghladata[i]<-ifelse(hlamerged$CWD[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="CWD", 
                                hlamerged$GFEs[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                                )], ifelse(hlamerged$CWD[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="NON-CWD", 
                                           hlamerged$GFEs[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1",  hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                                           )], ""))} 
    
    return(bigdawghladata)} else {print("Error: Unrecognized filename suffix. Stopping BDgenotypeconversion()")}}

hladata<-HLA_data
convertedlist<-sapply(colnames(atlas[,2:9]),function(x) NULL)

for(i in 1:length(mergedlist)){
convertedlist[[i]]<-BDgenotypeconversion(hladata, "/Users/liviatran/Desktop/ltmasterscoding/Allelelist.3310.txt", mergedlist[[i]])}


###########END SCRIPT  



bigdawghladata<-HLA_data


hlamerged<-cwdID("/Users/liviatran/Desktop/ltmasterscoding/Allelelist.3310.txt", mergedlist[[1]])
hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")
hlalist<-list()            
for (i in 3:ncol(bigdawghladata)){
  hlalist[[i]]<-strsplit(bigdawghladata[,i],":")} 
for(i in 3:ncol(bigdawghladata)){
  bigdawghladata[i]<-paste(sapply(hlalist[[i]], "[", 1), sapply(hlalist[[i]], "[", 2), sep=":")
  bigdawghladata[,i]<-paste(colnames(bigdawghladata[i]),bigdawghladata[,i],sep="*")}
for (i in 3:ncol(bigdawghladata)){
  bigdawghladata[i]<-ifelse(hlamerged$CWD[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="CWD", 
                            hlamerged$GFEs[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                            )], ifelse(hlamerged$CWD[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="NON-CWD", 
                                       hlamerged$GFEs[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1",  hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                                       )], ""))} 



