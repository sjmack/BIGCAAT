setwd("~/Desktop/ltmasterscoding")

#for file merging
multiFileread <- function(filepath, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  hlafiles<-lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  hlafiles<-Reduce(hlafiles, head(hlafiles, nrow(hlafiles)-(clip*length(filenames))), hlafiles)
  return(hlafiles)}

multiFileread("/Users/liviatran/Desktop/ltmasterscoding/HLA", skip=3, clip=1)
View(x[[1]])





hladf<-read.csv("hlaAdownload.csv", header=FALSE, stringsAsFactors=FALSE, skip=3, col.names= c("allelename", "gfes"))
hladf<-hladf[-seq(nrow(hladf),nrow(hladf)),]

#splits hladf with each respective gene feature getting its own column
hlasplit<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hladf[,2],"w")), nrow=nrow(hladf), byrow=T), stringsAsFactors=FALSE)) [2,],"-")
splitdf<-data.frame(matrix(unlist(hlasplit), nrow=nrow(hladf), byrow=T),stringsAsFactors=FALSE)

colnames(splitdf)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                      "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR")





