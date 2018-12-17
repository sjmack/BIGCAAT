#Functions for Variant Amino Acid Extraction from Locus Specific Exons
###This script is only for HLA-A -- future adaptations will be made for other HLA loci
#V 0.1
#By: L Tran
#12/13/18

##tested on A, B, and C

##Source Files

#Allelelist.3310.txt
#obtained from https://github.com/ANHIG/IMGTHLA/blob/Latest/Allelelist.3310.txt
#Various versions of Allelelist.___.txt are available depending on the version of HLA alleles used
#This file is a list of documented HLA alleles with their allele IDs
#needed for the alleleListfile parameter in cwdID

#MS_Eur.txt
#Patient HLA haplotypes in BIGDAWG format


########

##REQUIRED PACKAGES 
require(data.table)
require(stringr)
require(BIGDAWG)


###FUNCTIONS:
#function to count spaces in between regions of interest
#to determine where  start for the alignment sequence 
countSpaces <- function(x){
  counter <- 0
  coll <- numeric()
  vec <- strsplit(x," ")[[1]]
  for(i in 1:length(vec)){
    if (vec[i]==""){
      counter <- counter+1
    }
    else{
      if (counter!=0) coll <- c(coll,counter)
      counter <- 1
    }
  }
  coll
}


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

#reads in genotype data  
gdata<-read.table("MS_EUR.txt", sep="\t", header=T, check.names = F, stringsAsFactors = F)

#sets blank cells to NA 
#if cells do not contain NA, locus names are pasted to the allele in the MS_file
for (i in 3:ncol(gdata)){
  gdata[gdata==""]<-NA
  gdata[[i]]<-ifelse(is.na(gdata[[i]])==FALSE, paste(colnames(gdata[i]),gdata[,i],sep="*"), NA)}

#creates empty variables for future for loops
start<-end<-alignment<-list()

#creates empty variables where each element is named after the used loci 

#empty variables for correspondence table 
merp<-corr_table<-cols<-downloaded_segments<-w<-alignment_positions<-alignment_length<-alignment_start<-prot_extractions<-refblock_number<-end_char<-space_diff<-
  
  #empty variables for exon_extractor function   
  geno_alleles<-AA_segments<-AA_aligned <-refexon<-pepsplit<-alignment<-exonlist<- sapply(loci, function(x) NULL)

loci="A"



for(i in 1:length(loci)){
  #downloads relevant locus alignment file -- readLines allows for space preservation, which is important in
  #finding where the alignment sequence starts 
  alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(loci[[i]],"_prot.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
  #alters alignment file by cutting out non-pertinent information 
  alignment[[loci[i]]] <- head(alignment[[i]],-3)
  alignment[[loci[i]]] <- tail(alignment[[i]],-7)
  
  #see countSpaces function at the end of this script
  #Counts difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing to find where
  #the alignment sequence actually starts 
  if(loci[[i]]=="A"){
    space_diff[[loci[i]]]<-((countSpaces(alignment[[i]][1])[2])+(countSpaces(alignment[[i]][3])[2]))-(countSpaces(alignment[[i]][2])[1])}
  if(loci[[i]]=="B"){
    space_diff[[loci[i]]]<-((countSpaces(alignment[[i]][1])[2])+(countSpaces(alignment[[i]][3])[2]))-(countSpaces(alignment[[i]][2])[1])}
  if(loci[[i]]=="C"){
    space_diff[[loci[i]]]<-((countSpaces(alignment[[i]][1])[2])+(countSpaces(alignment[[i]][3])[2]))-(countSpaces(alignment[[i]][2])[1])}
  
  #if the loci is DRB1, space_diff is 0
  if(loci[[i]]=="DRB"){
    space_diff[[loci[[i]]]]<-((countSpaces(alignment[[i]][1])[2])+(countSpaces(alignment[[i]][3])[2]))-(countSpaces(alignment[[i]][2])[1])
  }
  
  
  #reduces repeated whitespace in alignment file and removes rows with empty values for proper
  #start and stop subsetting 
  alignment[[loci[i]]] <-str_squish(alignment[[i]])
  alignment[[loci[i]]] <-alignment[[i]][-which(alignment[[i]] == "")]
  
  #determines positions of "Prot" and the end of that reference block segment
  start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[i]]))
  
  end[[loci[i]]] <- as.numeric(c(start[[i]][2:length(start[[i]])]-1,length(alignment[[i]])))
  
  
  #counts number of characters in the very last allele to add onto the last Prot enumeration block
  #to obtain end length 
  end_char[[loci[i]]]<-nchar(sapply(strsplit(gsub(" ", "", sub(" ", "~", str_squish(tail(alignment[[i]], 1)))), "~"), "[", 2))-1
  
  
  #extracts rows with "Prot" and reference sequence position information 
  #extracts only relevant reference sequence positions
  #NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
  #as the actual sequence start will always be 1 
  for (j in 1:length(start[[i]])){
    
    prot_extractions[[i]][j]<-strsplit(alignment[[i]][start[[i]][j]], " ")
    
    refblock_number[[i]][j]<-as.numeric(sapply(prot_extractions[[i]][j], "[", 2))}
  
  
  #determines the alignment start by adding -30 to the difference between white spaces found above 
  alignment_start[[loci[i]]]<-refblock_number[[i]][1]+space_diff[[i]]
  
  
  #determines alignment length based on the last Prot enumeration + the number of characters in the last row 
  alignment_length[[loci[i]]]<-as.numeric(tail(refblock_number[[i]], 1))+end_char[[i]]
  
  #pastes alignment_start to alignment_length together in sequential order, with inDels accounted for 
  #captures output as "w"
  if(loci[i]=="A"){
    w$A<-capture.output(cat(alignment_start[[i]]:alignment_length[[i]]))}
  if(loci[i]=="B"){
    w$B<-capture.output(cat(alignment_start[[i]]:alignment_length[[i]]))}
  if(loci[[i]]=="C"){
    w$C<-capture.output(cat(alignment_start[[i]]:(as.numeric(tail(refblock_number[[i]],1))-(alignment_start[[i]])), paste("inDel", seq(1:6), sep=""), (as.numeric(tail(refblock_number[[i]],1))-alignment_start[[i]]+1):alignment_length[[i]]))}
  if(loci[[i]]=="DRB"){
    w$DRB<-capture.output(cat(alignment_start[[i]]:alignment_length[[i]]))}
  
  
  ##stopped here -- figured out DRB discrepancy between corr table and alignment position 
  #splits string formed by cat for separate character variables
  alignment_positions[[i]]<-as.character(unlist(strsplit(w[[i]], " ")))
  
  #eliminates "0", as the alignment sequence from ANHIG does not contain 0
  alignment_positions[[i]]<-alignment_positions[[i]][-which(alignment_positions[[i]] == 0)]
  
  
  #closes all white space in the alignment file, except for the white space separating the allele and peptide sequence
  alignment[[loci[i]]] <-paste(substr(alignment[[i]],1,regexpr(" ",text = alignment[[i]],fixed = TRUE)), gsub(" ","",substr(alignment[[i]],regexpr(" ",text = alignment[[i]],fixed = TRUE),nchar(alignment[[i]]))),sep = "")
  
  #string splits at white spaces to yield allele and peptide sequences
  alignment[[loci[i]]]  <- strsplit(alignment[[i]]," ", fixed=T)
  
  #binds the previously split strings by row and renames columns to "alleles" and "pepseq"
  alignment[[loci[i]]] <- do.call(rbind,alignment[[i]])
  
  alignment[[i]][which(alignment[[i]][,1]==alignment[[i]][,2]),2] <- ""
  
  colnames(alignment[[i]])<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")
  
  #due to ANHIG formatting, cases where an allele contains newly reference peptide sequences will not 
  #contain the same number of rows as previous reference peptide blocks
  #this for loop is invoked to add "."for all other alleles for each character in the newly reference peptide
  #to preserve structural integrity 
  for(k in 1:length(start[[i]])){
    if(nrow(alignment[[i]][start[[i]][k]:end[[i]][k],])!=nrow(alignment[[i]][start[[i]][1]:end[[i]][1],])){
      x<-as.data.frame(alignment[[i]][,1][start[[i]][1]:end[[i]][1]][-c(1,2)], stringsAsFactors = F)
      colnames(x)<-paste(loci[[i]], "alleles", sep="_")
      x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(alignment[[i]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
      y<-data.frame(tail(alignment[[i]],1), stringsAsFactors = F)
      x$pepseq[match(y[,1], x[,1])]<-y$pepseq
      alignment[[i]]<-as.matrix(rbind(head(alignment[[i]], -1), x))
      start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[i]]))
      end[[loci[i]]] <- as.numeric(c(start[[i]][2:length(start[[i]])]-1,nrow(alignment[[i]])))}}
  
  #if a locus has extra formatting, resulting in unqeual rows, start and end will be updated to reflect subsetting
  #if a locus has no extra formatting, start and end will remain the same, as procured by earlier code
  for(e in 1:length(start[[i]])){
    AA_segments[[i]]<-cbind(AA_segments[[i]], alignment[[i]][start[[i]][e]:end[[i]][e],])}
  
  #removes first two rows containing AA position and "Prot"
  AA_segments[[i]] <- AA_segments[[i]][-c(1,2),]
  
  #designates columns to be combined as every other so allele names are not included
  #in pasting all the amino acid sequences together 
  cols<-seq(0, ncol(AA_segments[[i]]), by=2)
  AA_segments[[i]]<-cbind(AA_segments[[i]][,1], apply(AA_segments[[i]][,cols], 1 ,paste, collapse = ""))
  
  #creates a new matrix with the number of columns equal to the number of characters in the reference sequence 
  corr_table[[i]]<-matrix(, nrow = 2, ncol = as.numeric(nchar(AA_segments[[i]][,2][1])))
  
  #contains actual sequence information
  corr_table[[i]][1,]<-(1:as.numeric(nchar(AA_segments[[i]][,2][1])))
  
  #contains alignment sequence information 
  corr_table[[i]][2,]<-alignment_positions[[i]]
  
  
  #string splits to extract locus in the allele name
  #assigns to new variable "AA_aligned"
  AA_aligned[[i]] <- as.matrix(do.call(rbind,strsplit(AA_segments[[i]][,1],"[*]")))
  
  #adds a new column of pasted locus and trimmed two field alleles to AA_aligned
  AA_aligned[[i]] <- cbind(AA_aligned[[i]], paste(AA_aligned[[i]][,1], apply(AA_aligned[[i]],MARGIN=c(1,2),FUN=GetField,Res=2)[,2], sep="*"))
  if(loci[[i]]=="DRB"){AA_aligned[[i]][,1]<-loci[[i]]}
  
  #binds AA_aligned and AA_segments -- renames columns 
  AA_segments[[i]] <- cbind(AA_aligned[[i]], AA_segments[[i]])
  colnames(AA_segments[[i]]) <- c("locus", "full_allele", "trimmed_allele", "allele_name", "AAsequence")
  
  #sets refexon to a reference peptide for each HLA locus based on the reference sequences in AA_segments 
  refexon[[i]] <- rbind(AA_segments[[i]][1,])[which(rbind(AA_segments[[i]][1,])[,"locus"]==loci[[i]]),'AAsequence']
  
  #splits AA_sequence column at every amino acid, resulting in a list of split amino acids for each row
  pepsplit[[i]] <- sapply(AA_segments[[i]][,"AAsequence"],strsplit,split="*")
  
  #fills in space with NA for alleles with premature termination to make it the same number of characters
  #as the reference sequence 
  pepsplit[[i]]<- lapply(pepsplit[[i]],function(x) c(x,rep("NA",nchar(refexon[[i]])-length(x))))
  
  #binds pep_split together by element in its previous list form by row
  pepsplit[[i]] <- do.call(rbind,pepsplit[[i]])
  
  #nullifies row names 
  rownames(pepsplit[[i]]) <- NULL
  
  #binds all columns together to form desired ouput, as described above
  AA_segments[[i]] <- cbind.data.frame(AA_segments[[i]][,1:4],pepsplit[[i]], stringsAsFactors=FALSE)
  
  #sets appropriate column names, where each amino acid has its respective position based on the
  #alignment sequence 
  colnames(AA_segments[[i]]) <- c("locus","allele","trimmed_allele","allele_name", corr_table[[1]][2,])
  
  #distributes  reference sequence from row 1
  #into all other rows, if they contain a "-"
  #amino acids with changes will not be impacted
  for(k in 5:ncol(AA_segments[[i]])) {
    AA_segments[[i]][,k][which(AA_segments[[i]][,k]=="-")] <- AA_segments[[i]][,k][1]}  
  
  #for loop for subsetting AA_sequences by matching exon start and end cells from AA_atlas
  #column names of AA_sequences, which are AA positions
  #subsets relevant amino acids, inputting them into a list
  #binds previous columns with locus, allele, trimmed allele, and allele name information
  
  #subsets first and last exons for each loci 
  exonlist[[i]][[1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][,5:match(as.numeric(AA_atlas[match(loci[[i]],names(AA_atlas))][[i]][[2]][[1]]-2), colnames(AA_segments[[i]]))])
  exonlist[[i]][[nrow(AA_atlas[[match(loci[[i]],names(AA_atlas))]])+1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][match(AA_atlas[[match(loci[[i]],names(AA_atlas))]][[2]][[length(AA_atlas[match(loci[[i]],names(AA_atlas))][[i]][[2]])]]:names(AA_segments[[i]][ncol(AA_segments[[i]])]), colnames(AA_segments[[i]]))])
  
  #for loop for subsetting N-1 exons
  for(j in 1:(nrow(AA_atlas[[match(loci[i],names(AA_atlas))]])-1)){
    exonlist[[i]][[j+1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][,match(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[j]], colnames(AA_segments[[i]])):match(as.numeric(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[j+1]]-1),colnames(AA_segments[[i]]))])}
  
  
  #for loop for subsetting exonlist alleles to only those found in genotype data
  #focuses on subsetting via the third column in exonlist, which consists of trimmed_allele data 
  for(d in 1:length(exonlist[[i]])){
    for(e in 1:2){
      
      #finds which exonlist alleles are present in genotype data alleles 
      if(loci[[i]]=="DRB")
      {geno_alleles[[e]]<-exonlist[[i]][[d]][,3][which(exonlist[[i]][[d]][,3] %in% gdata[which(colnames(gdata)%in%paste(loci[[i]], "1", sep="")==TRUE)][,e]==TRUE)]}
      
      else{
        geno_alleles[[e]]<-exonlist[[i]][[d]][,3][which(exonlist[[i]][[d]][,3] %in% gdata[which(colnames(gdata)%in%loci[[i]]==TRUE)][,e]==TRUE)]
      }}}
  
  #merges both sets of unique alleles found in exonlist and gets rid of duplicates 
  geno_alleles<-unique(append(geno_alleles[[1]], geno_alleles[[2]]))
  
  #creates a variable geno_exonlist, with the number of elements equal to how many exons there are for an allele
  geno_exonlist<-sapply(exonlist[[1]], function(x) NULL)
  
  #reads in 3310 HLA alleles 
  HLA_alleles<-read.csv("Allelelist.3310.txt", header=TRUE, stringsAsFactors = FALSE, skip=6,sep=",")
  
  #compiles a list of CWD alleles and inserts them into a new variable
  CWDalleles<-CWDverify()
  
  pastedAAseq<-cols<-all_gdata<-genotype_variants<-sapply(exonlist[[1]], function(x) NULL)
  nonCWD_checked<-singleAA_exon<-singleAA_alleles<-NULL
  
  #subsets exonlist alleles to those found in genotype data and inserts them into a new list
  #geno_exonlist
  for(d in 1:length(exonlist[[i]])){
    geno_exonlist[[d]]<-subset(exonlist[[i]][[d]], exonlist[[i]][[d]][,3]%in%geno_alleles)
    geno_exonlist[[d]]<-cbind.data.frame("accessions"=HLA_alleles$AlleleID[match(geno_exonlist[[d]]$allele_name, HLA_alleles$Allele)], geno_exonlist[[d]], stringsAsFactors=FALSE)
    geno_exonlist[[d]]<-cbind.data.frame("CWD"=ifelse(geno_exonlist[[d]]$accessions %in% CWDalleles$Accession, "CWD", "NON-CWD"), geno_exonlist[[d]], stringsAsFactors=FALSE)
    
    #subsets geno_exonlist to only containing CWD alleles via accession number
    #and stores it to a new variable, all_gdata
    #NOTE: all g_data will be a master copy of all variants of genotype data alleles
    if(any(geno_exonlist[[d]]$CWD=="CWD")){
      all_gdata[[d]]<-na.omit(geno_exonlist[[d]][geno_exonlist[[d]]$accessions%in%CWDalleles$Accession,])}
    
    #compares whether all truncated alleles in all_gdata are in geno_alleles
    #returns truncated alleles that are not CWD, but that are present in geno_alleles
    nonCWDtrunc<-cbind(geno_alleles%in%all_gdata[[d]]$trimmed_allele, geno_alleles)[which(cbind(geno_alleles, geno_alleles%in%all_gdata[[d]]$trimmed_allele)==FALSE)]
    
    #obtains non-CWD genotype variants in the genotype dataset
    for(b in 1:length(nonCWDtrunc)){
      genotype_variants[[d]][[b]]<-subset(geno_exonlist[[d]], geno_exonlist[[d]]$trimmed_allele==nonCWDtrunc[[b]])
      
      #if the non-CWD allele only has one variant, bind it to all_gdata
      if(nrow(genotype_variants[[d]][[b]])==1){all_gdata[[d]]<-rbind(all_gdata[[d]],genotype_variants[[d]][[b]])}
      
      #if the non-CWD allele has more than one variant, extract number of amino acid columns
      #present for a given exon 
      if(nrow(genotype_variants[[d]][[b]])>1){
        cols[[d]]<-7:length(genotype_variants[[d]][[b]])
        
        #if an exon for a non-CWD allele has more than one amino acid column, paste all the columns together to obtain
        #the amino acid sequence which is stored in pastedAAseq
        #pastedAAseq is evaluated to find which allele variant has the most complete sequence by counting the number of
        #character, omitting * (notation for unknown amino acid)
        #the allele with the most compelte sequence is bound to all_gdata
        if(length(cols[[d]])>1){
          pastedAAseq[[d]]<-apply(genotype_variants[[d]][[b]][ , cols[[d]]] , 1 , paste , collapse = "" )
          all_gdata[[d]]<-rbind(all_gdata[[d]], genotype_variants[[d]][[b]][names(pastedAAseq[[d]][which.max(nchar(gsub("[*^]","",pastedAAseq[[d]])))]),])}
        
        #if an exon for a non-CWD allele has one amino acid column (i.e. exon 8 for HLA-A), store it into a separate
        #variable, singleAA_alleles
        if(length(cols[[d]])==1){
          singleAA_exon[[b]]<-genotype_variants[[d]][[b]][ncol(genotype_variants[[d]][[b]])==7]
          singleAA_alleles<-singleAA_exon[lapply(singleAA_exon, length)>0]}}}
    
    #evaluates whether a variant amino acid is present and subsets it to nonCWD_checked if there is one
    #otherwise, if nonCWDchecked only contains *, use *
    for(c in 1:length(singleAA_alleles)){
      if(any(singleAA_alleles[[c]][7:length(singleAA_alleles[[c]])]!="*")==TRUE) {nonCWD_checked[[c]]<-subset(singleAA_alleles[[c]], singleAA_alleles[[c]][7:length(singleAA_alleles[[c]])]!="*")[1,]}
      if(any(singleAA_alleles[[c]][7:length(singleAA_alleles[[c]])]!="*")==FALSE){nonCWD_checked[[c]]<-subset(singleAA_alleles[[c]], singleAA_alleles[[c]][7:length(singleAA_alleles[[c]])]=="*")[1,]}
    }
    
    #binds narrowed down non-CWD alleles for one amino acid exons and inputs it back IF there is a one columned amino acid
    #if not, nothing happens 
    if(length(cols[[d]])==1)
    {all_gdata[[d]]<-rbind(all_gdata[[d]][ncol(all_gdata[[d]])==7], rbind(nonCWD_checked[[1]], nonCWD_checked[[2]]))}}
  
  
  #creates a new variable, position_parsed, with pre-defined elements based on
  #column names in AA_segments (i.e. position in the peptide sequence)
  position_parsed<-sapply(colnames(AA_segments[[i]][,5:ncol(AA_segments[[i]])]), function(x) NULL)
  
  #for loop to extract only variant amino acids and input them into their respective element positions
  #in position_parsed 
  #extracts only variant amino acids, discounting NA and unknown alleles (*)
  for(a in 1:length(all_gdata)){
    for(b in 1:length(7:ncol(all_gdata[[a]]))){
      position_parsed[match(colnames(all_gdata[[a]][7:ncol(all_gdata[[a]])]), names(position_parsed))][[b]]<-unique(subset(all_gdata[[a]][c(5,b+6)], (all_gdata[[a]][b+6]!=all_gdata[[a]][,b+6][1]) & (all_gdata[[a]][b+6] != "*") & (all_gdata[[a]][b+6] != "NA")))}}
  
}



position_parsed[sapply(position_parsed, nrow)>0]