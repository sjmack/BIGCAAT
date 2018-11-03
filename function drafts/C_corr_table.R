###Correspondence table for HLA-C
#By: Livia Tran
#11/2/2018
#v 0.2


loci="C"

#requires stringr package 
require(stringr)

#creates empty variables based on "loci" to be later filled in by for loops 
end<-corr_table<-downloaded_segments<-alignment_positions<-alignment_length<-alignment_start<-downloaded<-prot_extractions<-diff<-prot_start<-end_char<-refblock_number<-sapply(loci, function(x) NULL)


for(i in 1:length(loci)) {
  
  #downloads relevant locus alignment file -- readLines allows for space preservation, which is important in
  #finding where the alignment sequence starts 
  downloaded[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(loci,"_prot.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
  
  #alters downloaded file by cutting out non-pertinent information 
  downloaded[[loci[i]]] <- head(downloaded[[i]],-3)
  downloaded[[loci[i]]] <- tail(downloaded[[i]],-7)
  
  #see countSpaces function at the end of this script
  #Counts difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing to find where
  #the alignment sequence actually starts 
  diff[[loci[i]]]<-(countSpaces(downloaded[[i]][2])[1])-(countSpaces(downloaded[[i]][1])[2])+1
  
  #reduces repeated whitespace in downloaded file and removes rows with empty values for proper
  #start and stop subsetting 
  downloaded[[loci[i]]] <-str_squish(downloaded[[i]])
  downloaded[[loci[i]]] <-downloaded[[i]][-which(downloaded[[i]] == "")]
  
  #determines positions of "Prot" and the end of that reference block segment
  prot_start[[loci[i]]]<-as.numeric(grep("Prot", downloaded[[i]]))
  end[[loci[i]]] <- as.numeric(c(prot_start[[i]][2:length(prot_start[[i]])]-1,length(downloaded[[i]])))
  
  #counts number of characters in the very last allele to add onto the last Prot enumeration block
  #to obtain end length 
  end_char[[loci[i]]]<-nchar(sapply(strsplit(gsub(" ", "", sub(" ", "~", str_squish(tail(downloaded[[i]], 1)))), "~"), "[", 2))-1
  
  #extracts rows with "Prot" and reference sequence position information 
  #extracts only relevant reference sequence positions
  #NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
  #as the actual sequence start will always be 1 
  for (j in 1:length(prot_start)){
    prot_extractions[[j]]<-downloaded[[i]][prot_start[[j]]]
    prot_extractions[[j]]<-strsplit(prot_extractions[[j]], " ")
    refblock_number[[j]]<-as.numeric(sapply(prot_extractions[[j]], "[", 2))
    
    #determines the alignment start by adding -30 to the difference between white spaces found above 
    alignment_start[[loci[i]]]<-refblock_number[[j]][1]+diff[[i]]
    
    #determines alignment length based on the last Prot enumeration + the number of characters in the last row 
    alignment_length[[loci[i]]]<-as.numeric(tail(refblock_number[[j]], 1))+end_char[[i]]
    
    #pastes alignment_start to alignment_length together in sequential order, with inDels accounted for 
    #captures output as "w"
    w<-capture.output(cat(alignment_start[[i]]:(as.numeric(refblock_number[[j]][4])-(alignment_start[[i]])), paste("inDel", seq(1:6), sep=""), (as.numeric(refblock_number[[i]][4])-alignment_start[[i]]+1):alignment_length[[i]]))
    
    #splits string formed by cat for separate character variables
    alignment_positions[[i]]<-as.character(unlist(strsplit(w, " ")))
    
    #eliminates "0", as the alignment sequence from ANHIG does not contain 0
    alignment_positions[[i]]<-alignment_positions[[i]][-which(alignment_positions[[i]] == 0)]}
  
  #closes all white space in the downloaded file, except for the white space separating the allele and peptide sequence
  downloaded[[loci[i]]] <-paste(substr(downloaded[[i]],1,regexpr(" ",text = downloaded[[i]],fixed = TRUE)), gsub(" ","",substr(downloaded[[i]],regexpr(" ",text = downloaded[[i]],fixed = TRUE),nchar(downloaded[[i]]))),sep = "")
  
  #string splits at white spaces to yield allele and peptide sequences
  downloaded[[loci[i]]]  <- strsplit(downloaded[[i]]," ", fixed=T)
  
  #binds the previously split strings by row and renames columns to "alleles" and "pepseq"
  downloaded[[loci[i]]] <- do.call(rbind,downloaded[[i]])
  colnames(downloaded[[i]])<-c(paste(loci[[i]], "alleles", sep=""), "pepseq")
  
  ##if the downloaded alignment file contains updates where alleles have extra amino acids, the row will not be equal
  #to the rows that previous start:end subsets contain
  #if all rows are not equal, thus containing new alleles with new amino acids, 
  #create a new block segment, where each allele contains "." for the number of characters the new allele contains
  #if all rows are equal, nothing will change 
  for(k in 1:length(prot_start[[i]])){
    if(nrow(downloaded[[i]][prot_start[[i]][k]:end[[i]][k],])!=nrow(downloaded[[i]][prot_start[[i]][1]:end[[i]][1],]))
    {x<-as.data.frame(downloaded[[i]][,1][prot_start[[i]][1]:end[[i]][1]][-c(1,2)], stringsAsFactors = F)
    colnames(x)<-paste(loci[[i]], "alleles", sep="")
    x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(downloaded[[i]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
    y<-data.frame(tail(downloaded[[i]],1), stringsAsFactors = F)
    x$pepseq[match(y[,1], x[,1])]<-y$pepseq
    downloaded[[i]]<-as.matrix(rbind(head(downloaded[[i]], -1), x))
    prot_start[[loci[i]]]<-as.numeric(which(downloaded[[i]][,1]=="Prot"))
    end[[loci[i]]] <- as.numeric(c(prot_start[[i]][2:length(prot_start[[i]])]-1,nrow(downloaded[[i]])))}}}

for(i in 1:length(loci)){
  for(k in 1:length(prot_start[[i]])){
    
  #subsets downloaded file by start and stop positions -- binds all subsets together, where the number of characters
  #in the reference sequence may be obtained to input as the end position for the correspondence table 
  downloaded_segments[[i]]<-cbind(downloaded_segments[[i]], downloaded[[i]][prot_start[[i]][k]:end[[i]][k],])}
  cols<-seq(0, ncol(downloaded_segments[[i]]), by=2)
  downloaded_segments[[i]]<-cbind(downloaded_segments[[i]][,1], apply(downloaded_segments[[i]][,cols], 1 ,paste, collapse = ""))
  rownames(downloaded_segments[[i]])=NULL
  
  #creates a new matrix with the number of columns equal to the number of characters in the reference sequence 
  corr_table[[i]]<-matrix(, nrow = 2, ncol = as.numeric(nchar(downloaded_segments[[i]][,2][3])))
  
  #contains actual sequence information
  corr_table[[i]][1,]<-(1:as.numeric(nchar(downloaded_segments[[i]][,2][3])))
  
  #contains alignment sequence information 
  corr_table[[i]][2,]<-alignment_positions[[i]]
}


##function to count spaces in between regions of interest
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


