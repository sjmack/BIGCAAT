#Extraction of core and non-core exon amino acid sequences 
#V 0.2
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
AA_segments <- AA_segments[-c(1,2),]

AA_segments<-cbind(AA_segments[,1],apply(AA_segments[,2:ncol(AA_segments)],MARGIN=1,paste,collapse=""))


