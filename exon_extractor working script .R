#Extraction of core and non-core exon amino acid sequences 
#V 0.1
#By: Livia Tran

#fetches "_prot.txt files from IMGTHLA for a specified HLA locus
filefetcher <- function(loci) {
  #downloads *_prot.txt alignment files

  # Get Locus Based Alignments
  for(i in 1:length(loci)) {
    locus <- loci[i]
    filename <- paste(locus,"_prot.txt",sep="")
    URL <- paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",FileName,sep="")
    download.file(URL,destfile = filename,method="libcurl")}}
  
#fetches _prot.txt alignment files -- files are input to current working directory 
GetFiles(c("A", "B", "C", "DPA1", "DPB1", "DRB1", "DRB3", "DRB4", "DRB5"))  

#a function to extract core exon and non-core exon amino acid sequences
exon_extractor <- function(locus){
name <- paste(locus,"_prot.txt",sep="")
alignment <- read.table(name,fill=T,header=F,stringsAsFactors=F,strip.white=T,colClasses="character")
alignment <- as.matrix(alignment[-c(1:2),]) #Remove Header
alignment <- as.matrix(alignment[-nrow(alignment),]) #Remove Footer
} 

x<-exon_extractor("A")  
View(x[- grep("Prot", x),])

