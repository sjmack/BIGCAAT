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
    dataset <- read.table(dataset,header = TRUE,sep="\t",check.names = TRUE)
  } else { colnames(dataset) <- make.names(colnames(dataset),unique = TRUE)}
  
  # Everything is stored in the strataSet list 
  stratSet <- data.frame(rep(NA,nrow(dataset)))
  
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
