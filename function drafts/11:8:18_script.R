#THIS IS A PENDING NEW VERSION OF EXON_EXTRACTOR BASED ON A NEW ATLAS BEING TESTED ON HLA-A
#REVISIONS WILL BE APPENDED TO THE THE CURRENT EXON_EXTRACTOR_WORKING SCRIPT UPON SJ MACK'S INPUT






##AA_atlas for HLA-A
AA_atlas<-list()
loci<-"A"

i=1:7
AA_atlas[[loci[1]]] <- data.frame(exon=paste(i, i+1, sep=":"), bound=c(1, 91, 183, 275, 314, 325, 341),stringsAsFactors = FALSE)


########

#required packages 
require(data.table)
require(stringr)
require(BIGDAWG)

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
corr_table<-cols<-downloaded_segments<-w<-alignment_positions<-alignment_length<-alignment_start<-prot_extractions<-refblock_number<-end_char<-space_diff<-
  
  #empty variables for exon_extractor function   
  geno_alleles<-AA_segments<-AA_aligned <-refexon<-pepsplit<-alignment<-exonlist<- sapply(loci, function(x) NULL)


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
  space_diff[[loci[i]]]<-(countSpaces(alignment[[i]][2])[1])-(countSpaces(alignment[[i]][1])[2])+1
  
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
    w$A<-capture.output(cat(alignment_start[[i]]:(as.numeric(refblock_number[[i]][4])-(alignment_start[[i]])), (as.numeric(refblock_number[[i]][4])-alignment_start[[i]]+1):alignment_length[[i]]))}
  if(loci[i]=="C"){
    if(loci[[i]]=="C"){
      w$C<-capture.output(cat(alignment_start[[i]]:(as.numeric(refblock_number[[i]][4])-(alignment_start[[i]])), paste("inDel", seq(1:6), sep=""), (as.numeric(refblock_number[[i]][4])-alignment_start[[i]]+1):alignment_length[[i]]))}    }
  
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
      alignment[[i]]<-as.matrix(rbind(head(alignment[[i]], -1), x))}
    
    
    #if/when all start and end position segments have the same number of rows,
    #starts and ends for each set of amino acid sequences are parsed out where the
    #output is a matrix with allele names and each subsequent set of amino acid sequences
    else{AA_segments[[i]]<-cbind(AA_segments[[i]], alignment[[i]][start[[i]][k]:end[[i]][k],])}}
  
  
  #removes first two rows containing AA position and "Prot"
  AA_segments[[i]] <- AA_segments[[i]][-c(1,2),]
  
  #designates columns to be combined as every other so allele names are not included
  #in pasting all the amino acid sequences together 
  cols<-seq(0, ncol(AA_segments[[i]]), by=2)
  AA_segments[[i]]<-cbind(AA_segments[[i]][,1], apply(AA_segments[[i]][,cols], 1 ,paste, collapse = ""))
  
  
  #creates a new matrix with the number of columns equal to the number of characters in the reference sequence 
  corr_table[[i]]<-matrix(, nrow = 2, ncol = as.numeric(nchar(AA_segments[[1]][,2][1])))
  
  #contains actual sequence information
  corr_table[[i]][1,]<-(1:as.numeric(nchar(AA_segments[[1]][,2][1])))
  
  #contains alignment sequence information 
  corr_table[[i]][2,]<-alignment_positions[[i]]
  
  #string splits to extract locus in the allele name
  #assigns to new variable "AA_aligned"
  AA_aligned[[i]] <- as.matrix(do.call(rbind,strsplit(AA_segments[[i]][,1],"[*]")))
  
  #adds a new column of pasted locus and trimmed two field alleles to AA_aligned
  AA_aligned[[i]] <- cbind(AA_aligned[[i]], paste(AA_aligned[[i]][,1], apply(AA_aligned[[i]],MARGIN=c(1,2),FUN=GetField,Res=2)[,2], sep="*"))
  
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
  exonlist[[i]][[1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][,5:match(as.numeric(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[1]]-2), colnames(AA_segments[[i]]))])
  exonlist[[i]][[length(AA_atlas[[i]][[1]])+1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][match(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[7]], colnames(AA_segments[[i]]))])
  
  
  #for loop for subsetting N-1 exons
  for(j in 1:(nrow(AA_atlas[[match(loci[i],names(AA_atlas))]])-1)){
    exonlist[[i]][[j+1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][,match(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[j]], colnames(AA_segments[[i]])):match(as.numeric(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[j+1]]-1),colnames(AA_segments[[1]]))])}
  
  #for loop for subsetting exonlist alleles to only those found in genotype data
  #focuses on subsetting via the third column in exonlist, which consists of trimmed_allele data 
  for(d in 1:length(exonlist[[i]])){
    for(e in 1:length(gdata[which(colnames(gdata)%in%exonlist[[i]][[d]][,1]==TRUE)])){
      
      #finds which exonlist alleles are present in genotype data alleles 
      geno_alleles[[e]]<-exonlist[[i]][[d]][,3][which(exonlist[[i]][[d]][,3] %in% gdata[which(colnames(gdata)%in%exonlist[[i]][[d]][,1]==TRUE)][,e]==TRUE)]}}
  
  #merges both sets of unique alleles found in exonlist and gets rid of duplicates 
  geno_alleles<-unique(append(geno_alleles[[1]], geno_alleles[[2]]))
  
  
  geno_exonlist<-sapply(exonlist[[1]], function(x) NULL)
  
  #subsets exonlist alleles to those found in genotype data and inserts them into a new list
  #geno_exonlist
  for(d in 1:length(exonlist[[i]])){
    for(e in 1:length(gdata[which(colnames(gdata)%in%exonlist[[i]][[d]][,1]==TRUE)])){
      geno_exonlist[[d]]<-subset(exonlist[[i]][[d]], exonlist[[i]][[d]][,3]%in%geno_alleles)}}
  
  #creates variables with the number of elements being equal to the length of the "actual" sequence
  position_parsed<-sapply(colnames(AA_segments[[i]][,5:ncol(AA_segments[[i]])]), function(x) NULL)
  
  #for loop to extract only variant amino acids and input them into their respective element positions
  #in position_parsed 
  #extracts only variant amino acids, discounting NA and unknown alleles (*)
  for(a in 1:length(geno_exonlist)){
    for(b in 1:length(5:ncol(geno_exonlist[[a]]))){
      position_parsed[match(colnames(geno_exonlist[[a]][5:ncol(geno_exonlist[[a]])]), names(position_parsed))][[b]]<-subset(geno_exonlist[[a]][c(1,3,b+4)], (geno_exonlist[[a]][b+4]!=geno_exonlist[[a]][,b+4][1]) & (geno_exonlist[[a]][b+4] != "*") & (geno_exonlist[[a]][b+4] != "NA"))}}
  
  #removes elements with no variant amino acids (i.e. a row length of 0)
  position_parsed<-position_parsed[sapply(position_parsed, nrow)>0]
  
}


#creates a new list with 2 columns, since each locus in gdata has 2 columns
variant_match<-sapply(1:2, function(x) NULL)

#further adds on to the previously created variable, where the final product is 
#a list with 2 elements
#where each element contains a list of 831 elements, i.e. the number of rows in gdata
#where each of the 831 elements contains a list of 145 elements, i.e. the length of position_parsed
for(i in 1:length(variant_match)){
  variant_match[[i]]<-sapply(nrow(gdata), function(x) NULL)
  for(j in 1:nrow(gdata)){
    variant_match[[i]][[j]]<-sapply(position_parsed, function(x) NULL)
  }
}

#only subsets data with where the locus column of position_parsed matches the column name of gdata
#subsets every single cell in corresponding matching columns to return variant amino acids
for(f in 1:length(position_parsed)){
  for(g in 1:length(gdata[which(colnames(gdata)%in%position_parsed[[f]][,1])])){
    for(h in 1:nrow(gdata)){
      variant_match[[g]][[h]][[f]]<-subset(position_parsed[[f]][3], position_parsed[[f]][2]==gdata[which(colnames(gdata)%in%position_parsed[[f]][,1])][g][h,])}
  }}

#removes elements with zero rows (i.e. no matches)
for(f in 1:length(position_parsed)){
  for(g in 1:length(gdata[which(colnames(gdata)%in%position_parsed[[f]][,1])])){
    for(h in 1:nrow(gdata)){
      variant_match[[g]][[h]]<-variant_match[[g]][[h]][sapply(variant_match[[g]][[h]], nrow)>0]}}}


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



