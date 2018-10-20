#exon_extractor
#V 0.5 - 10/19/18
#By: Livia Tran

###NOTE: this script is in its early stages and is using only the A_prot.txt and B_prot.txt files 
#for code execution

#requires data.table package
require(data.table)

#downloads appropriate loci_prot.txt files from ANHIG and reads them in
#each individual locus will be an element alignment, an empty list created to
#store protein sequences for each locus 
#skips first 5 lines containing citation information
#clips last three lines containing blank spaces and impertinent information
#removes white space between peptide sequences for all loci
#preserves white space between allele name and peptide sequence for strsplit later on
#as.matrix dimensions for strsplit later on 
loci=c("A", "B")
alignment<-list()
for(i in 1:length(loci)) {
  alignment[[loci[i]]] <- fread(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/3330/alignments/",paste(loci[i],"_prot.txt",sep=""),sep=""), fill=T, sep="\t", skip=7, header=F,strip.white=T,colClasses="character")
  alignment[[loci[i]]] <- as.matrix(head(alignment[[i]],-3)) 
  alignment[[loci[i]]] <-as.matrix(paste(substr(alignment[[i]],1,regexpr(" ",text = alignment[[i]],fixed = TRUE)), gsub(" ","",substr(alignment[[i]],regexpr(" ",text = alignment[[i]],fixed = TRUE),nchar(alignment[[i]]))),sep = ""))}


#creates empty variables for the next for loop to input start/end values for each HLA locus
start<-list()
end<-list()

#splits allele names from peptide sequences for each locus
#binds  previously split list rowwise -- output is an allele name column and the AA sequences column
#for each locus
#finds columns where the allele name is equal to the AA sequences - occurs when expression is affected
#due to some mutation
#defines start as where "Prot" is present
#since Prot is at the start of the next set of AA sequences
#defines end as the second element in start -1 row to signify the end of that segment
#includes the final row of the matrix
for(i in 1:length(alignment)){
alignment[[loci[i]]]  <- strsplit(alignment[[i]]," ", fixed=T)
alignment[[loci[i]]]  <- as.matrix(do.call(rbind,alignment[[i]]))
alignment[[i]][which(alignment[[i]][,1]==alignment[[i]][,2]),2] <- ""
start[[loci[i]]]<-as.numeric(which(alignment[[i]][,1]=="Prot"))
end[[loci[i]]] <- as.numeric(c(start[[i]][2:length(start[[i]])]-1,nrow(alignment[[i]])))
}

#for loop for parsing out starts and ends for each set of amino acid sequences
#output is a matrix with allele names and each subsequent set of amino acid sequences
AA_segments <- sapply(names(alignment),function(x) NULL)

for(i in 1:length(alignment)){
  for(j in 1:length(start[[i]])){
    AA_segments[[i]]<-cbind(AA_segments[[i]], alignment[[i]][start[[i]][j]:end[[i]][j],])}}

#for loop for removing headers, which contains AA position and "Prot"
#designates columns to be combined as every other so allele names are not included
#in pasting all the amino acid sequences together 
for(i in 1:length(AA_segments)){
AA_segments[[i]] <- AA_segments[[i]][-c(1,2),]
cols<-seq(0, ncol(AA_segments[[i]]), by=2)
AA_segments[[i]]<-cbind(AA_segments[[i]][,1], apply(AA_segments[[i]][,cols], 1 ,paste, collapse = ""))}

#for loop for string split to extract locus in the allele name
#assigns to new variable "AA_aligned"
#trims allele name down to two fields -- based on GetField function from BIGDAWG general functions
#combines previous AA_segments matrix with AA_aligned matrix, for a total of 5 columns
#AA_aligned extracts locus information to input to AA_segments 
#assigns column names based on what each column contains 
AA_aligned<-sapply(names(alignment),function(x) NULL)

for(i in 1:length(AA_segments)){
AA_aligned[[i]] <- do.call(rbind,strsplit(AA_segments[[i]][,1],"[*]"))
AA_segments[[i]] <- cbind(cbind(AA_aligned[[i]],apply(AA_aligned[[i]],MARGIN=c(1,2),FUN=GetField,Res=2)[,2]),AA_segments[[i]])
colnames(AA_segments[[i]]) <- c("locus","full_allele","trimmed_allele","allele_name","AAsequence")}

#sets refexon to a reference peptide for each HLA locus based on the reference sequences in
#AA_segments 
refexon<- sapply(names(alignment),function(x) NULL)
for(i in 1:length(AA_segments)){
refexon[[i]] <- rbind(AA_segments[[i]][1,])[which(rbind(AA_segments[[i]][1,])[,"locus"]==loci[[i]]),'AAsequence']}

#sets start as 1 for each loci
exon_start<-sapply(names(alignment),function(x) NULL)

for(i in 1:length(AA_segments)){
exon_start[[i]]<-1}

#splits AA_sequence column at every amino acid, resulting in a list of split amino acids for each row
#binds pep_split together by element in its previous list form by row
#nullifies row names 
#binds the first 4 columns of AA_segments with amino acid positions from pepsplit 
pepsplit<-sapply(names(alignment),function(x) NULL)

for(i in 1:length(AA_segments)){
pepsplit[[i]] <- sapply(AA_segments[[i]][,"AAsequence"],strsplit,split="*")
pepsplit[[i]]<- lapply(pepsplit[[i]],function(x) c(x,rep("NA",nchar(refexon[[i]])-length(x))))
pepsplit[[i]] <- do.call(rbind,pepsplit[[i]])
rownames(pepsplit[[i]]) <- NULL
AA_segments[[i]] <- cbind(AA_segments[[i]][,1:4],pepsplit[[i]])
colnames(AA_segments[[i]]) <- c("locus","allele","trimmed_allele","allele_name",seq(exon_start[[i]],ncol(AA_segments[[i]][,5:ncol(AA_segments[[i]])])+exon_start[[i]]-1))
}

#for loop for distributing the reference sequence from row 1
#into all other rows, if the contain a "-"
#amino acids with changes will not be impacted
for(i in 1:length(AA_segments)){
  for(j in 5:ncol(AA_segments[[i]])) {
    AA_segments[[i]][,j][which(AA_segments[[i]][,j]=="-")] <- AA_segments[[i]][,j][1]
  }}

#adds in *00:00 alleles to account for alleleic absences for each locus previously defined by loci
#at the beginning of the script 
#turns matrix into a dataframe
for(i in 1:length(AA_segments)){
AA_segments[[i]] <- data.frame(rbind(c(loci[[i]],"00:00:00:00","00:00",paste(loci[[i]],"*00:00:00:00",sep=""), rep(".",ncol(AA_segments[[i]])-4)),
AA_segments[[i]]), stringsAsFactors = FALSE, check.names = F)}

#AA_atlas.R is used as a guide for determining start and stop positions for each exon for a given HLA locus
load("AA_atlas.rda")

#for loop for subsetting AA_sequences by matching exon start and end cells from AA_atlas
#column names of AA_sequences, which are AA positions
#subsets relevant amino acids, inputting them into a list
#binds previous columns with locus, allele, trimmed allele, and allele name information 
exonlist<-sapply(names(alignment),function(x) NULL)
for(i in 1:length(AA_segments)){
  for(j in 1:nrow(AA_atlas[[i]]))
exonlist[[i]][[j]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][match(AA_atlas[[i]][j,2],colnames(AA_segments[[i]])):match(AA_atlas[[i]][j,3],colnames(AA_segments[[i]]))])}

View(exonlist[[1]][[1]])

