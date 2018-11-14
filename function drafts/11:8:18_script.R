#THIS IS A PENDING NEW VERSION OF EXON_EXTRACTOR BASED ON A NEW ATLAS BEING TESTED ON HLA-A
#REVISIONS WILL BE APPENDED TO THE THE CURRENT EXON_EXTRACTOR_WORKING SCRIPT UPON SJ MACK'S INPUT


##AA_atlas for HLA-A
AA_atlas<-list()
loci<-"A"

i=1:7
AA_atlas[[loci[1]]] <- data.frame(exon=c(paste("exon", paste(i, i+1, sep=":"), sep="_")), stringsAsFactors = FALSE)
AA_atlas[[loci[[1]]]]$bound <-list(1, 91, 183, 275, 314, 325, 341)



########
AA_segments<-AA_aligned <-refexon<-pepsplit<-alignment<-exonlist<- sapply(loci, function(x) NULL)

#creates empty variables for the next for loop to input start/end values for each HLA locus
start<-end<-list()

###BEGIN FOR LOOP 

#downloads appropriate loci_prot.txt files from ANHIG and reads them in
#each individual locus will be an element alignment, an empty list created to
#store protein sequences for each locus 
#skips first 5 lines containing citation information
#clips last three lines containing blank spaces and impertinent information
#removes white space between peptide sequences for all loci
#preserves white space between allele name and peptide sequence for strsplit later on
#as.matrix dimensions for strsplit later on 
#splits allele names from peptide sequences for each locus
#binds  previously split list rowwise -- output is an allele name column and the AA sequences column
#for each locus
#finds columns where the allele name is equal to the AA sequences - occurs when expression is affected
#due to some mutation
#defines start as where "Prot" is present
#since Prot is at the start of the next set of AA sequences
#defines end as the second element in start -1 row to signify the end of that segment
#includes the final row of the matrix
alignment<-list()
for(i in 1:length(loci)) {
  alignment[[loci[i]]] <- fread(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(loci[i],"_prot.txt",sep=""),sep=""), fill=T, sep="\t", skip=7, header=F,strip.white=T,colClasses="character")
  alignment[[loci[i]]] <- as.matrix(head(alignment[[i]],-3)) 
  alignment[[loci[i]]] <-as.matrix(paste(substr(alignment[[i]],1,regexpr(" ",text = alignment[[i]],fixed = TRUE)), gsub(" ","",substr(alignment[[i]],regexpr(" ",text = alignment[[i]],fixed = TRUE),nchar(alignment[[i]]))),sep = ""))
  alignment[[loci[i]]]  <- strsplit(alignment[[i]]," ", fixed=T)
  alignment[[loci[i]]]  <- as.matrix(do.call(rbind,alignment[[i]]))
  alignment[[i]][which(alignment[[i]][,1]==alignment[[i]][,2]),2] <- ""
  colnames(alignment[[i]])<-c(paste(loci[[i]], "alleles", sep=""), "pepseq")
  start[[loci[i]]]<-as.numeric(which(alignment[[i]][,1]=="Prot"))
  end[[loci[i]]] <- as.numeric(c(start[[i]][2:length(start[[i]])]-1,nrow(alignment[[i]])))
  for(j in 1:length(start[[i]])){
    if(nrow(alignment[[i]][start[[i]][j]:end[[i]][j],])!=nrow(alignment[[i]][start[[i]][1]:end[[i]][1],]))
    {x<-as.data.frame(alignment[[i]][,1][start[[i]][1]:end[[i]][1]][-c(1,2)], stringsAsFactors = F)
    colnames(x)<-paste(loci[[i]], "alleles", sep="")
    x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(alignment[[i]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
    y<-data.frame(tail(alignment[[i]],1), stringsAsFactors = F)
    x$pepseq[match(y[,1], x[,1])]<-y$pepseq
    alignment[[i]]<-as.matrix(rbind(head(alignment[[i]], -1), x))
    start[[loci[i]]]<-as.numeric(which(alignment[[i]][,1]=="Prot"))
    end[[loci[i]]] <- as.numeric(c(start[[i]][2:length(start[[i]])]-1,nrow(alignment[[i]])))}}}

#for loop for parsing out starts and ends for each set of amino acid sequences
#output is a matrix with allele names and each subsequent set of amino acid sequences
for(i in 1:length(loci)){
  for(j in 1:length(start[[i]])){
    AA_segments[[i]]<-cbind(AA_segments[[i]], alignment[[i]][start[[i]][j]:end[[i]][j],])}}

#for loop for removing headers, which contains AA position and "Prot"
#designates columns to be combined as every other so allele names are not included
#in pasting all the amino acid sequences together 
#string splits to extract locus in the allele name
#assigns to new variable "AA_aligned"
#trims allele name down to two fields -- based on GetField function from BIGDAWG general functions
#combines previous AA_segments matrix with AA_aligned matrix, for a total of 5 columns
#AA_aligned extracts locus information to input to AA_segments 
#assigns column names based on what each column contains 
#sets refexon to a reference peptide for each HLA locus based on the reference sequences in
#AA_segments 
#splits AA_sequence column at every amino acid, resulting in a list of split amino acids for each row
#binds pep_split together by element in its previous list form by row
#nullifies row names 
#binds the first 4 columns of AA_segments with amino acid positions from pepsplit 
for(i in 1:length(loci)){
  AA_segments[[i]] <- AA_segments[[i]][-c(1,2),]
  cols<-seq(0, ncol(AA_segments[[i]]), by=2)
  AA_segments[[i]]<-cbind(AA_segments[[i]][,1], apply(AA_segments[[i]][,cols], 1 ,paste, collapse = ""))
  AA_aligned[[i]] <- do.call(rbind,strsplit(AA_segments[[i]][,1],"[*]"))
  AA_segments[[i]] <- cbind(cbind(AA_aligned[[i]],apply(AA_aligned[[i]],MARGIN=c(1,2),FUN=GetField,Res=2)[,2]),AA_segments[[i]])
  colnames(AA_segments[[i]]) <- c("locus","full_allele","trimmed_allele","allele_name","AAsequence")
  refexon[[i]] <- rbind(AA_segments[[i]][1,])[which(rbind(AA_segments[[i]][1,])[,"locus"]==loci[[i]]),'AAsequence']
  pepsplit[[i]] <- sapply(AA_segments[[i]][,"AAsequence"],strsplit,split="*")
  pepsplit[[i]]<- lapply(pepsplit[[i]],function(x) c(x,rep("NA",nchar(refexon[[i]])-length(x))))
  pepsplit[[i]] <- do.call(rbind,pepsplit[[i]])
  rownames(pepsplit[[i]]) <- NULL
  AA_segments[[i]] <- cbind.data.frame(AA_segments[[i]][,1:4],pepsplit[[i]], stringsAsFactors=FALSE)
  alignment_seq<-seq(as.numeric(corr_table[[1]][2,][1]), as.numeric(corr_table[[1]][2,][1])+ncol(AA_segments[[i]][,5:ncol(AA_segments[[i]])])+1-1)
  alignment_seq<-alignment_seq[-which(alignment_seq == 0)]
  colnames(AA_segments[[i]]) <- c("locus","allele","trimmed_allele","allele_name", alignment_seq)}

#for loop for distributing the reference sequence from row 1
#into all other rows, if the contain a "-"
#amino acids with changes will not be impacted
for(i in 1:length(loci)){
  for(j in 5:ncol(AA_segments[[i]])) {
    AA_segments[[i]][,j][which(AA_segments[[i]][,j]=="-")] <- AA_segments[[i]][,j][1]
  }}



#for loop for subsetting AA_sequences by matching exon start and end cells from AA_atlas
#column names of AA_sequences, which are AA positions
#subsets relevant amino acids, inputting them into a list
#binds previous columns with locus, allele, trimmed allele, and allele name information 

for(i in 1:length(AA_segments)){
  for(j in 1:(nrow(AA_atlas[[match(loci[i],names(AA_atlas))]])-1)){
    exonlist[[i]][[1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][,5:match(as.numeric(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[1]]-2), colnames(AA_segments[[i]]))])
    exonlist[[i]][[j+1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[i]][,match(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[j]], colnames(AA_segments[[i]])):match(as.numeric(AA_atlas[match(loci[i],names(AA_atlas))][[i]][[2]][[j+1]]-1),colnames(AA_segments[[i]]))])
    exonlist[[i]][[length(AA_atlas[[1]][[1]])+1]]<-cbind(AA_segments[[i]][,1:4], AA_segments[[1]][match(AA_atlas[match(loci[i],names(AA_atlas))][[1]][[2]][[7]], colnames(AA_segments[[1]]))])}}


#creates variables with the number of elements being equal to the length of the "actual" sequence
position_parsed<-split_alleles<-sapply(colnames(AA_segments[[1]][,5:ncol(AA_segments[[1]])]), function(x) NULL)

#for loop to extract only variant amino acids and input them into their respective element positions
#in position_parsed 
for(k in 1:length(exonlist[[1]])){
for(i in 1:length(position_parsed[match(colnames(exonlist[[1]][[k]][5:ncol(exonlist[[1]][[k]])]), names(position_parsed))])){
  for(j in (5:ncol(exonlist[[1]][[k]]))[i]){
    position_parsed[match(colnames(exonlist[[1]][[k]][5:ncol(exonlist[[1]][[k]])]), names(position_parsed))][[i]]<-subset(exonlist[[1]][[k]][c(4,j)], exonlist[[1]][[k]][j]!=exonlist[[1]][[k]][,j][1])}}}

#reads in MS_file 
MS_file<-read.table("MS_EUR.txt", sep="\t", header=T, check.names = F, stringsAsFactors = F)

#sets blank cells to NA 
#if cells do not contain NA, locus names are pasted to the allele in the MS_file
for (i in 3:ncol(MS_file)){
MS_file[MS_file==""]<-NA
MS_file[[i]]<-ifelse(is.na(MS_file[[i]])==FALSE, paste(colnames(MS_file[i]),MS_file[,i],sep="*"), NA)}

#creates a variable with the number of elements being equal to the length of the "actual" sequence
#creates another list under each position with 6 more elements, where each element is a locus found in 
#MS_files 
parsed_list<-lapply(split_alleles, function(x) append(x, list(A=NA, C=NA, B=NA, DRB1=NA, DQB1=NA, DPB1=NA)))

#trims full allele names to a resolution of 2 
for(i in 1:length(parsed_list)){
split_alleles[[i]]<-strsplit(sapply(strsplit(position_parsed[[i]][,1], "\\*"), "[", 2), ":")
position_parsed[[i]][,1]<-paste(sapply(strsplit(position_parsed[[i]][,1], "\\*"), "[", 1), paste(sapply(split_alleles[[i]], "[", 1), sapply(split_alleles[[i]], "[", 2), sep=":"), sep="*")}

is.data.frame(position_parsed[[1]][1])

##for loop to sort through position_parsed data
##where alleles found in MS_files are extracted
#those alleles are converted to their respective variant amino acid
#and input into their appropriate locus identifier within the amino acid position element 
for(i in 1:length(parsed_list)){
    for(j in 1:length(parsed_list[[i]])){
      for(k in 3:length(na.omit(match(colnames(MS_file), names(parsed_list[[i]][j]))))){
    ifelse(names(parsed_list[[i]][j])%in%colnames(MS_file[k]), parsed_list[[i]][which(names(parsed_list[[i]][j])%in%colnames(MS_file[k])==TRUE)]<-cbind(position_parsed[[i]][,1][match(MS_file[,k], position_parsed[[i]][,1])], position_parsed[[i]][,2][match(MS_file[,k], position_parsed[[i]][,1])]), NA)}}} 



View(position_parsed[[1]][,2][match(MS_file[,3], position_parsed[[1]][,1]),drop=FALSE])

