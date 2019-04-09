###Sequence Feature Isolation 
#Module 1 Functions 
##By: Livia Tran 
#V 1.15
#September 30, 2018

##This script aims to isolate pre-determined feature groups for a given locus in its GFE notation.

#HLA data may be compared to zeroed-out GFEs to obtain appropriate
#GFE notations, depending on the feature sequence of interest. 
#Zeroed-out GFE notations refer to all fields other than the feature(s) of interest to contain zeros.
#Analysis may then be carried out by BIGDAWG to better understand distribution of categories and 
#its effects on phenotypes (case or control).
#Through two main functions, individual feature sequences or feature groups can be isolated in obtaining
#GFE information for BIGDAWG formatted data.
#customGFEgenerator() generates a set of custom HLA GFE notations based on a user-defined ‘atlas’ that 
#directs the ‘zeroing-out’ of specific features and feature groups.
#BIGDAWG_GFEanalyzer() converts each allele in a BIGDAWG formatted genotype dataset into its corresponding 
#GFE notation, as defined by a given set of custom GFE notations, followed by analysis in BIGDAWG
#for outputs of analyzed locus data.

#Individual features and feature group isolations are based off a user-made atlas. The atlas is based off
#a framework, which is also user-made.
#The framework is a list of HLA loci, with the names of individual gene features derived from
#GFE notations (i.e. 5'UTR, Exon 1, Intron 1...)
#The atlas relies on the framework for locus information, which is put in the first column of the atlas.
  #The atlas is meant to be a "road-map" for users to define what feature group or individual features
  #they want compare. This atlas consists of individual features previously defined in the framework
  #as well as feature groups informed by structural domains described in
  #LeFranc et. al, 2005, "IMGT unique numbering for MHC groove G-Domain and MHC superfamily (MhcSF) G-LIKE-DOMAIN".
  #Feature groups that contain one exon (i.e. an individual feature) 
  #are combined name-wise (see exon 1 and Leader peptide, exon 8 and C-domain). Coordinates for a desired
  #feature of interest allows all other gene features to be zeroed-out, thus focusing on the feature of interest.
  #Each column in the atlas is a "map" that defines the location of each feature or feature group for a 
  #given locus.
  
##Source-Files:
#BSG files
  #Generated from Bx12 Search Generator, which compiles data from the Neo4J database. BSG files
  #are downloadable by locus, where each file consists of an allele name for a given locus and its
  #GFE notations. 11 total HLA files where downloaded from Bx12 Search Generator in .csv format.
  #These are used in multiFileread().
#Allelelist.3310.txt
  #obtained from https://github.com/ANHIG/IMGTHLA/blob/Latest/Allelelist.3310.txt
  #Various versions of Allelelist.___.txt are available depending on the version of HLA alleles used
  #This file is a list of documented HLA alleles with their allele IDs
  #needed for the alleleListfile parameter in cwdID

###loads pre-made RDA files 
###framework.rda is a list consisting of feature identification for each field of a GFE
###atlas2.0.rda is a dataframe that gives the position for selected features in a GFE as well as individual features
load("/Users/liviatran/ltmasterscoding/framework.rda")
load("/Users/liviatran/ltmasterscoding/atlas2.0.rda")

#stringr package is necessary for str_replace_all function
require(stringr)

#BIGDAWG package is necessary for sample BIGDAWG formatted HLA data, and 
#BIGDAWG function is used later for Locus analysis 
require(BIGDAWG)

############FUNCTIONS 
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


#for file merging
filemerge <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<-Reduce(function(x,y) {merge(x,y,all=TRUE)}, datalist)
  datalist<-head(datalist,nrow(datalist)-(clip*length(filenames)))
  return(datalist)}

#for identifying CWD alleles and identifying each allele's unique allele ID 
##optional usage for alleles with CWD alleles (i.e. HLA)
#users looking at other genes that do not have CWD documentation need not use this function
cwdID <- function(allelelistFile,gfeDirPath) {
  hlaacc<-read.csv(allelelistFile, header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")
  hladf<-filemerge(gfeDirPath, c("allelename", "gfe"), skip=3, clip=1)
  hladf$alleleID<-hlaacc$AlleleID[match(hladf$allelename, paste("HLA",hlaacc$Allele,sep="-"))] 
  cwdalleles<-CWDverify()
  hladf$CWD<-ifelse(hladf$alleleID %in% cwdalleles$Accession, "CWD", "NON-CWD")
  return(hladf)}

##function for reading all BSG files in, combines them into a single, multi-locus dataframe
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
  if(any(unlist(atlas[x,y])==0) == FALSE) {dataframe[,unlist(atlas[x,y])] <- GFEsplitDF[unlist(atlas[x,y])]}
  return(dataframe)}

#function for pasting all GFE fields together
#adds an additional row to all loci for each sequence feature
#inputs a *00:00 allele, with a maxed out GFE for each field 
GFEpaster<-function(seqfeaturelist, seqfeatureposition, BSGdf, BSGdfposition){
  cols<- 1:length(seqfeaturelist[[seqfeatureposition]])
  GFEs<-str_replace_all(paste(gsub("[*].*", "", unlist(BSGdf[[BSGdfposition]][[1]])),paste("w", apply(seqfeaturelist[[seqfeatureposition]][,cols] , 1 , paste, collapse="-"))), " ", "")
  featureselectDF<-cbind.data.frame(BSGdf[[BSGdfposition]][[1]], GFEs, stringsAsFactors=FALSE)
  colnames(featureselectDF)[1] <- colnames(BSGdf[[BSGdfposition]][,])[1]
  return(featureselectDF)}

#function to convert BIGDAWG formatted data into its GFE component 
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
    for (i in 3:ncol(bigdawghladata)){
      hlalist[[i]]<-strsplit(bigdawghladata[,i],":")} 
    for(i in 3:ncol(bigdawghladata)){
      bigdawghladata[i]<-ifelse(is.na(hlalist[[i]])==FALSE, paste(sapply(hlalist[[i]], "[", 1), sapply(hlalist[[i]], "[", 2), sep=":"), NA)    
      bigdawghladata[i]<-ifelse(is.na(bigdawghladata[[i]])==FALSE, paste(colnames(bigdawghladata[i]),bigdawghladata[,i],sep="*"), NA)}
    for (i in 3:ncol(bigdawghladata)){
      bigdawghladata[i]<- ifelse(is.na(bigdawghladata[[i]])==FALSE, 
                                 ifelse(hlamerged$CWD[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"),)]=="CWD",
                                        hlamerged$GFEs[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))], ifelse(hlamerged$CWD[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="NON-CWD", hlamerged$GFEs[match(bigdawghladata[,i], paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))], NA)),NA)}
    return(bigdawghladata)} else {print("Error: Unrecognized filename suffix. Stopping BDgenotypeconversion()")}}


###########BEGIN SCRIPT for customGFEgenerator()

####Function for generating a custom list of HLA-GFE tables, stored in mergedlist at the end of the function 
customGFEgenerator<-function(GFEfiledirectory, columnnames, skip, clip){
hladf<-multiFileread(GFEfiledirectory, columnnames, skip=skip, clip=clip)

#splits GFE notations -- first at "w", then at "-"
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w", fixed = TRUE)), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-", fixed=TRUE))

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
seqfeatureslist <- sapply(colnames(atlas[,2:length(atlas)]),function(x) NULL)

#for loop for inputting sequence feature isolations for all loci in list format 
for (k in 1:length(seqfeatureslist)) {
  for (i in 1:nrow(atlas)){
    seqfeatureslist[[k]][[atlas$locus[i]]]<-featureselect(i, k+1, splitDFs[[i]])
    seqfeatureslist[[k]][[atlas$locus[i]]][nrow(seqfeatureslist[[k]][[atlas$locus[i]]])+1,] <- as.data.frame(rbind(rep(0,length(seqfeatureslist[[k]][[atlas$locus[i]]]))))
    if(any(unlist(atlas[i,k+1])==0) == FALSE) {seqfeatureslist[[k]][[atlas$locus[i]]][nrow(seqfeatureslist[[k]][[atlas$locus[i]]]),][,unlist(atlas[i,k+1])]<-unlist(rep(max(as.numeric(as.character(unlist(splitDFs[[i]]))))+1,length(framework[[i]])))[unlist(atlas[i,k+1])]}
      }}

#appends new row to BSG files to account for *00:00 alleles 
for (i in 1:nrow(atlas)){
  hladf[[i]][nrow(hladf[[i]]) + 1,] = list(as.character(str_replace_all(paste(atlas$locus[i], "*00:00"), " ", "")), NA)}

#makes an empty list of 8 total features, names of elements defined by atlas column names
isolatedlist <- sapply(colnames(atlas[,2:length(atlas)]),function(x) NULL)

#for loop for pasting all sequence feature fields together
#inputs into isolatedlist
for(i in 1:nrow(atlas)){
  for(k in 1:length(seqfeatureslist)){
    isolatedlist[[k]][[atlas$locus[i]]]<-GFEpaster(seqfeatureslist[[k]], i, hladf, i)}}

##makes an empty list of 8 total features, names of elements defined by atlas column names
mergedlist<-sapply(colnames(atlas[,2:length(atlas)]),function(x) NULL)

#recursively merges all dataframes for each sequence feature together from isolated list
#inputs into merged list 
for(i in 1:length(isolatedlist)){
  mergedlist[[i]]<-as.data.frame(Reduce(function(x,y) {merge(x,y,all=TRUE)}, isolatedlist[[i]]))}

return(mergedlist)
}

###########END SCRIPT for customGFEgenerator()


###########BEGIN SCRIPT for BIGDAWG_GFEanalyzer()
##function to convert BIGDAWG formatted data into its GFE counterpart based on feature group or sequence feature desired

#parameters to note:
#info is a logical parameter that allows the user to view map options (i.e. to look at all possible options for
#feature groups or individual sequence features)
#info is defaulted to FALSE, where the user must input the desired feature sequence or feature group into the mapname argument
#if a user specifies info=T, a message with map options is ouput to the console
#the mapname parameter identifies the desired map of the atlas, which directs the function to return
#GFE notations for that desired map
#from the custom merged data set from customGFEgenerator(), and outputs GFE notations for that specific
#mapname
#if the user defines mapname=="all", all feature sequences and feature groups present in the atlas
#are sequentially used to define the conversion of HLA alleles to atlas-specified custom GFE 
#notations in a BIGDAWG-formatted genotype dataset.
#if the user does not specify anything for mapname, the default is always "all"
#return is a logical parameter as an option for users to view their custom
#feature group or individual feature isolated GFE notations before they are
#analyzed by BIGDAWG 
BIGDAWG_GFEanalyzer<-function(BIGDAWGgenotypedata, alleleListfiles, mergedcustomdata=custom_mergeddata, mapname="all", info=F, return=F){
  #logical parameter for info
  if(info==TRUE){cat(paste("The following ‘maps’ are available in the atlas:",paste(colnames(atlas)[2:ncol(atlas)],collapse=" "),sep="\n"))}
  else{
    #default start and end to 0 
    start<-end<-0
     
    #creates an empty list if only one individual feature/feature group is desired
    convertedlist <- list()
    
    #captures the passed BIGDAWGgenotypedata
    datasetName <- match.call()[2]
    
    #sets class to character 
     class(datasetName) <- "character"
     
##Logical Parameters
     #if map name = all, start=1, end is the length of names of merged custom data, proceed=T
     #else
     #if mapname is in names of mergedcustomdata, start=mapname, end=mapname, proceed=T
     #if map name is not in names of mergedcustomdata or all, print an error message and do NOT proceed
    if(mapname=="all"){
      start=1
      end=length(names(mergedcustomdata))
      proceed<-TRUE}
     
     else{
       if(mapname%in%names(mergedcustomdata)==TRUE){
       start<-end<-match(mapname,names(mergedcustomdata))
       proceed<-TRUE}
     
       else{
         cat("Invalid map name. Please set info=T to view map names, or input 'all' to use all map names.", sep="\n")
         proceed<-FALSE}}
   ###end of logical parameters 
     
  ###if proceed--   
     #for loop for converting BIGDAWG formatted data into its GFE equivalent, where start and end
     #are defined by the logical parameters above, depending on what the user chooses for mapname
     #convertedlist is named after the input BIGDAWGgenotype data - name of map name 
     if(proceed){ 
     for (i in start:end){
      convertedlist[[paste(datasetName,names(mergedcustomdata)[i],sep="-")]]<-BDgenotypeconversion(BIGDAWGgenotypedata, alleleListfiles, mergedcustomdata[[i]])
       cat(paste("*** Analyzing",names(mergedcustomdata[i]),"dataset ***","\n",sep=" "))
       BIGDAWG(convertedlist[[paste(datasetName,names(mergedcustomdata)[i],sep="-")]], HLA=F, Run.Tests ="L", Verbose=FALSE)}
    
      #if return =T, the user may view the converted GFE equivalent of the BIGDAWG formatted data 
     if(return){
        return(convertedlist)}
     }}}

############END SCRIPT for BIGDAWG_GFEanalyzer()

#### BDStrat -- simple allele stratification of BIGDAWG formatted datasets  
####            Steven J. Mack September 13 - 18, 2018 v0.1

###             BDStrat() accepts a BIGDAWG formatted genotype data file and generates a list containing a pair of subset 
###             dataframes (named locus*allelles-positive and locus*alleles-negative) in which one dataframe contains all 
###             subjects with a specified allele, and the second dataframe contains all subjects without the specified 
###             allele. These dataframes can then be passed to BIGDAWG for analyses stratified on the specified allele. 

###             Current version (v0.1) is limited to stratifying on multiple alleles at a single locus
###             Current version assumes the BIGDAWG dataset is tab-delimited

# Examples 
# stratified <- BDStrat("MS_EUR.txt","DRB1",c("15:01","11:04"))
# stratified2 <- BDStrat(HLA_data,"DRB1",c("15:01:01:01","11:04:01"))

BDStrat <- function(dataset,locus,alleles){
  
  #for(i in length(locus)){ ## hook for a future (more complicated) version that stratifies on multiple loci
  locus <- c(locus,paste(locus,"1",sep="."))
  #}
  
  # Since all loci are duplicated, check.names = TRUE generates "locus","Locus.1" name pairs for read files
  # and make.names(x,unique=TRUE) does the same for data frames
  if(!is.data.frame(dataset)) { 
    dataset <- read.table(dataset,header = TRUE,sep="\t",check.names = TRUE, stringsAsFactors = FALSE)
  } else { colnames(dataset) <- make.names(colnames(dataset),unique = TRUE)}
  
  # Everything is stored in the strataSet list 
  stratSet <- data.frame(rep(NA,nrow(dataset)), stringsAsFactors = FALSE)
  
  # Identify the rows containing each target allele for each locus column
  for(i in 1:length(alleles)){
    stratSet <- cbind(stratSet,dataset[[locus[1]]]==alleles[i],dataset[[locus[2]]]==alleles[i])
  }
  
  # Identify the rows containing any target allele
  for(i in 1:nrow(stratSet)){
    stratSet[i,1] <- any(unlist(stratSet[i,2:((2*length(alleles))+1)]))
  }
  
  # Split the parent dataset into two stratified subsets
  posStrat <- dataset[stratSet[,1]==TRUE,]
  negStrat <- dataset[stratSet[,1]==FALSE,]
  
  # Add them as elements of the stratPair list, named with the selected or excluded alleles
  stratPair <- list()
  stratPair[[paste(locus[1],"*",paste(unlist(alleles),collapse=paste("+",locus,"*",sep="")),"-positive",sep="")]] <- posStrat
  stratPair[[paste(locus[1],"*",paste(unlist(alleles),collapse=paste("+",locus,"*",sep="")),"-negative",sep="")]] <- negStrat
  
  # Strip suffixes from column names
  colnames(stratPair[[1]]) <- gsub(".1","",colnames(stratPair[[1]]),fixed=TRUE)
  colnames(stratPair[[2]]) <- gsub(".1","",colnames(stratPair[[2]]),fixed=TRUE)
  
  #return object
  stratPair
}

####EXAMPLES

##tests out customGFEgenerator
#"/Users/liviatran/Desktop/ltmasterscoding/HLA" is a list of BSG files for all HLA loci 
custom_mergeddata<-customGFEgenerator("/Users/liviatran/Desktop/ltmasterscoding/HLA", columnnames = c("allelename", "gfe"), skip=3, clip=1)

#save custom_mergeddata -- only generate again of BSG files are updated 
save(custom_mergeddata,file="custom_mergeddata.rda")

#tests out BIGDAWGGFEanalyzer
BIGDAWG_GFEanalyzer("MS_EUR.txt","Allelelist.3310.txt", mapname = "threeUTR")

#stratifies HLA_data to negatively and positively associated MS alleles
stratified<-(BDStrat("MS_EUR.txt","DRB1","15:01"))

#for loop for converting both negatively and positively associated MS alleles to their
#GFE notations followed by analysis by BIGDAWG
for(i in 1:length(stratified)){
cat(paste("*** Analyzing",names(stratified[i]),"dataset ***","\n",sep=""))
BIGDAWG_GFEanalyzer(stratified[[i]], "Allelelist.3310.txt", mapname="threeUTR")}



