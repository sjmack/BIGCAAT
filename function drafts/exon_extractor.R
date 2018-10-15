#exon_extractor
#V 0.4
#By: Livia Tran


###NOTE: this script is in its early stages and is using only the A_prot.txt file for code execution
##areas where the locus is hard-coded will be fixed later on when the code is adapted for
#universal usage for all HLA loci 
#fetches "_prot.txt files from IMGTHLA for a specified HLA locus
#downloads directly to working directory
#additionally reads the .txt file in, storing protein information for each locus as an element
#in a list named "alignment" 
#strips white space
#fill=T for blank rows
loci="A"
alignment<-list()
for(i in 1:length(loci)) {
    locus <- loci[i]
    filename <- paste(locus,"_prot.txt",sep="")
    URL <- paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",filename,sep="")
    download.file(URL,destfile = filename, method="libcurl")
    name <- paste(loci,"_prot.txt",sep="")
    alignment[[name[i]]]<- read.table(name[i],fill=T,header=F,stringsAsFactors=F,strip.white=T,colClasses="character", sep="\t")
  }

#removes footer at the end of the text file
alignment <- as.matrix(alignment[[1]][-nrow(alignment[[1]]),]) 

#removes white space between peptide sequences for all elements in alignment 
for(i in 1:nrow(alignment)){
  alignment[i]<-paste(substr(alignment[i],1,regexpr(" ",text = alignment[i],fixed = TRUE)[1]),gsub(" ","",substr(alignment[i],regexpr(" ",text = alignment[i],fixed = TRUE)[1],nchar(alignment[i]))),sep="")}

#splits column 1 of alignment at the white space between
#allele and amino acid sequence to delineate allele name from AA sequences
alignment <- strsplit(alignment[,1]," ")

#binds the previously split list rowwise -- output is an allele name column and the AA sequences column
alignment <- as.matrix(do.call(rbind,alignment))

#finds columns where the allele name is equal to the AA sequences - occurs when expression is affected
#due to some mutation
alignment[which(alignment[,1]==alignment[,2]),2] <- ""

#defines start as where "Prot" is present
#since Prot is at the start of the next set of AA sequences
start<-which(alignment[,1]=="Prot")

#defines end as the second element in start -1 row to signify the end of that segment
#includes the final row of the matrix
end <- c(start[2:length(start)]-1,nrow(alignment))

#for loop for parsing out starts and ends for each set of AA sequences
#output is a matrix with allele names and each subsequent set of AA sequences 
AA_segments<-NULL
for(i in 1:length(start)){
  AA_segments <- cbind(AA_segments,alignment[start[i]:end[i],])}

#removes header, which contains AA position and "Prot"
AA_segments <- AA_segments[-c(1,2),]

#designates columns to be combined as every other so allele names are not included
#in pasting all the amino acid sequences together 
cols<-seq(0, ncol(AA_segments), by=2)
AA_segments<-cbind(AA_segments[,1], apply(AA_segments[,cols], 1 ,paste, collapse = ""))

#string split to extract locus in the allele name
#assigns to new variable
AA_aligned <- do.call(rbind,strsplit(AA_segments[,1],"[*]"))

#trims allele name down to two fields -- based on GetField function from BIGDAWG general functions
#combines previous AA_segments matrix with AA_aligned matrix, for a total of 5 columns
#assigns column names based on what each column contains 
AA_segments <- cbind(cbind(AA_aligned,apply(AA_aligned,MARGIN=c(1,2),FUN=GetField,Res=2)[,2]),AA_segments)
colnames(AA_segments) <- c("locus","full_allele","trimmed_allele","allele_name","AAsequence")

#commented out for further development when testing multiple loci
##for later usage when testing out on multiple alleles -- ensures locus specific rows
#AA_segments <- AA_segments[which(AA_segments[,'locus']=="A"),]

#creates a reference table data frame with the first row of AA_segments, which contains 
#reference information
#binds with a new column row named ref.start, which contains information for where 
#an exon starts 
reftable<-as.data.frame(cbind(rbind(AA_segments[1,]), ref.start=1), stringsAsFactors = FALSE)

####in this case, the locus is hard-coded to A, as this code is still in development before testing
#it out on all HLA loci 
#sets RefExon to a Reference Peptide, depending on locus desired 
RefExon <- reftable[which(reftable[,"locus"]=="A"),'AAsequence']

#sets RefStart to a locus specific start position 
RefStart <- as.numeric(reftable[which(reftable[,'locus']=="A"),'ref.start'])

#sets RefAllele to a locus specific reference allele 
RefAllele <- reftable[which(reftable[,'locus']=="A"),'allele_name']

#designates start position in the full reference peptide sequence in AA_segments
#by finding what position it starts at based on RefExon
exon_start <- as.numeric(regexpr(RefExon,AA_segments[1,'AAsequence']))

#designates end position in the full reference peptide sequence in AA_segments
#by adding exon start to the number of characters in RefExon -1 
exon_end <- (exon_start + nchar(RefExon))-1

#subsets reference peptide sequences in AA_segments based on found exon start and end in previous lines
AA_segments[,'AAsequence'] <- substr(AA_segments[,'AAsequence'],exon_start,exon_end)

#splits AA_sequence column at every amino acid, resulting in a list of split amino acids for each row
pep_split <- sapply(AA_segments[,'AAsequence'],strsplit,split="*")

#fills all allele rows with "-" until the length of the RefExon
#rows with less characters are likely due some mutation 
pep_split<- lapply(pep_split,function(x) c(x,rep("NA",nchar(RefExon)-length(x))))

#binds pep_split together by element in its previous list form by row
#nullifies row names 
pep_split <- do.call(rbind,pep_split)
rownames(pep_split) <- NULL

#binds the first 4 columns of AA_segments with amino acid positions from pep_split 
AAsequences <- cbind(AA_segments[,1:4],pep_split)

#renames each column
#with each AA getting its own column name, depending on its position in the original peptide sequence
colnames(AAsequences) <- c("locus","allele","trimmed_allele","allele_name",seq(RefStart,ncol(AAsequences[,5:ncol(AAsequences)])+RefStart-1))

#for loop for distributing the reference sequence from row 1
#into all other rows, if the contain a "-"
#amino acids with changes will not be impacted
for(i in 5:ncol(AAsequences)) {
  x <- AAsequences[,i]
  x[which(x=="-")] <- x[1]
  AAsequences[,i] <- x
}

#NOTE: locus is hard-coded here -- will be fixed later
#adds in *00:00 alleles to account for alleleic absences 
#turns matrix into a dataframe
AAsequences <- as.data.frame(rbind(c("A","00:00:00:00","00:00",paste("A","*00:00:00:00",sep=""), rep(".",ncol(AAsequences)-4)),
                                   AAsequences), stringsAsFactors = FALSE)
               

#AA_atlas.R is used as a guide for determining start and stop positions for each exon for a given HLA locus
load("AA_atlas.rda")

#creates empty list
exonlist<-list()

#for loop for subsetting AA_matrix by matching exon start and end cells from AA_atlas
#column names of AA_matrix, which are AA positions
#subsets relevant amino acids, inputting them into a list
#binds previous columns with locus, allele, trimmed allele, and allele name information 
for(i in 1:nrow(AA_atlas)){
exonlist[[AA_atlas$exon[i]]]<-cbind(AAsequences[,1:4], AAsequences[,match(AA_atlas[i,2],colnames(AAsequences)):match(AA_atlas[i,3],colnames(AAsequences))])}


View(exonlist[[1]])

