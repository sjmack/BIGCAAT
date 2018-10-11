#Extraction of core and non-core exon amino acid sequences 
#V 0.3
#By: Livia Tran

#fetches "_prot.txt files from IMGTHLA for a specified HLA locus
#downloads directly to working directory
loci=c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5")
  for(i in 1:length(loci)) {
    locus <- loci[i]
    filename <- paste(locus,"_prot.txt",sep="")
    URL <- paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",filename,sep="")
    download.file(URL,destfile = filename, method="libcurl")}

###NOTE: this script is in its early stages and is using only the A_prot.txt file for code execution
##areas where the locus is hard-coded will be fixed later on when the code is adapted for
#universal usage for all HLA loci 
name <- paste("A","_prot.txt",sep="")

#reads in loci_prot.txt files, with the locus defined by "name" 
#currently only looking at HLA-A
#strips white space
#fill=T for blank rows
alignment <- read.table(name,fill=T,header=F,stringsAsFactors=F,strip.white=T,colClasses="character", sep="\t")

#removes footer at the end of the text file
alignment <- as.matrix(alignment[-nrow(alignment),]) 

#replaces the first instance of white space with ~ for split of allele names from reference amino
#acid sequences
alignment[,1] <- sapply(alignment[,1],FUN=sub,pattern=" ",replacement="~")

#replaces all instances of white space with no white space 
alignment[,1] <- sapply(alignment[,1],FUN=gsub,pattern=" ",replacement="")

#splits column 1 of alignment at ~ to delineate allele name from AA sequences
alignment <- strsplit(alignment[,1],"~")

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
  if(increments <- end[i]-start[i])
  {AA_segments <- cbind(AA_segments,alignment[start[i]:end[i],])}}

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
AA_aligned <- cbind(AA_aligned,apply(AA_aligned,MARGIN=c(1,2),FUN=GetField,Res=2)[,2])

#combines previous AA_segments matrix with AA_aligned matrix, for a total of 5 columns
#assigns column names based on what each column contains 
AA_segments <- cbind(AA_aligned,AA_segments)
colnames(AA_segments) <- c("locus","full_allele","trimmed_allele","allele_name","AAsequence")

#inputs the full reference sequence in row 1 of AA_segments into a variable called
#reference_sequence 
reference_sequence <- AA_segments[1,]

##for later usage when testing out on multiple alleles -- ensures locus specific rows
AA_segments <- AA_segments[which(AA_segments[,'locus']=="A"),]

#binds AA_segments with reference sequence to input the reference sequence without impacting
#the locus specific allele containing the full reference sequence
#nullifies row names
AA_segments <- rbind(reference_sequence,AA_segments)
AA_segments[1,1:4]  <- "reference sequence"
rownames(AA_segments) <- NULL

#calls upon RefTab, a reference table from BIGDAWG
#contains Locus, Reference Locus, Reference Allele, Reference Exons, Reference Peptide,
#and Reference Start
RefTab <- BIGDAWG::ExonPtnList$RefExons[1,]

#changes Reference Exons from only core exons to include all exons
RefTab$Reference.Exons<-list(c(1, 2, 3, 4, 5, 6, 7, 8))

#changes Reference Peptide to include all amino acid sequences from Exon 1 - end of Exon 
RefTab$Reference.Peptide<-"MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAVHAAEQRRVYLEGRCVDGLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWELSSQPTIPIVGIIAGLVLLGAVITGAVVAAVMWRRKSSDRKGGSYTQAASSDSAQGSDVSLTACKV"

RefTab$Reference.Start<-"1"

####in this case, the locus is hard-coded to A, as this code is still in development before testing
#it out on all HLA loci 
#sets RefExon to a Reference Peptide, depending on locus desired 
RefExon <- RefTab[which(RefTab[,"Locus"]=="A"),'Reference.Peptide']

#sets RefStart to a locus specific start position 
RefStart <- as.numeric(RefTab[which(RefTab[,'Locus']=="A"),'Reference.Start'])

#sets RefAllele to a locus specific reference allele 
RefAllele <- RefTab[which(RefTab[,'Locus']=="A"),'Reference.Allele']

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
pep_split<- lapply(pep_split,function(x) c(x,rep("-",nchar(RefExon)-length(x))))

#binds pep_split together by element in its previous list form by row
#nullifies row names 
pep_split <- do.call(rbind,pep_split)
rownames(pep_split) <- NULL

#binds the first 4 columns of AA_segments with amino acid positions from pep_split 
AA_matrix <- cbind(AA_segments[,1:4],pep_split)

#renames each column
#with each AA getting its own column name, depending on its position in the original peptide sequence
colnames(AA_matrix) <- c("locus","allele","trimmed_allele","allele_name",paste("position",seq(RefStart,ncol(AA_matrix[,5:ncol(AA_matrix)])+RefStart-1),sep="."))

#for loop for distributing the reference sequence from row 1
#into all other rows, if the contain a "-"
#amino acids with changes will not be impacted
for(i in 5:ncol(AA_matrix)) {
  x <- AA_matrix[,i]
  x[which(x=="-")] <- x[1]
  AA_matrix[,i] <- x
}

#removes reference sequence 
AA_matrix <- AA_matrix[-1,]

#NOTE: locus is hard-coded here -- will be fixed later
#adds in *00:00 alleles to account for alleleic absences 
AA_matrix <- rbind(c("A","00:00:00:00","00:00",paste("A","*00:00:00:00",sep=""),rep("^",ncol(AA_matrix)-4)),
                   AA_matrix)

#creates empty list
exonlist<-list()

#AA_atlas.R is used as a guide for determining start and stop positions for each exon for a given HLA locus
load("AA_atlas.rda")

#for loop for subsetting AA_matrix by matching exon start and end cells from AA_atlas
#column names of AA_matrix, which are AA positions
#subsets relevant amino acids, inputting them into a list
for(i in 1:nrow(AA_atlas)){
exonlist[[AA_atlas$exon[i]]]<-AA_matrix[2:nrow(AA_matrix),match(AA_atlas[i,2],colnames(AA_matrix)):match(AA_atlas[i,3],colnames(AA_matrix))]}

View(exonlist[[1]])

