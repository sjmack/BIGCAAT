setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-C

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

hlaC<-data.frame(matrix(unlist(hlasplit[3]), byrow=T, nrow=length(hlasplit[[3]])), stringsAsFactors = FALSE)
colnames(hlaC)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                  "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR")


template<-hlaC
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,17)))

##UTRS
C_UTR_template<-template
C_UTR_template[,unlist(atlas[3,2])] <- hlaC[unlist(atlas[3,2])]

#CORE EXONS
C_CE_template<-template
C_CE_template[,unlist(atlas[3,3])] <- hlaC[unlist(atlas[3,3])]

#CYTO EXONS
C_CYTO_template<-template
C_CYTO_template[,unlist(atlas[3,4])] <- hlaC[unlist(atlas[3,4])]

#LEADER PEPTIDE
C_LPEP_template<-template
C_LPEP_template[,unlist(atlas[3,5])] <- hlaC[unlist(atlas[3,5])]

#C-DOMAIN
C_CDOM_template<-template
C_CDOM_template[,unlist(atlas[3,6])] <- hlaC[unlist(atlas[3,6])]

#INTRONS
C_INTRONS_template<-template
C_INTRONS_template[,unlist(atlas[3,7])] <- hlaC[unlist(atlas[3,7])]

#TMEXONS
C_TMEXONS_template<-template
C_TMEXONS_template[,unlist(atlas[3,8])] <- hlaC[unlist(atlas[3,8])]


#NON-CORE
C_NONCORE_template<-template
C_NONCORE_template[,unlist(atlas[3,9])] <- hlaC[unlist(atlas[3,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}