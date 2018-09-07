
setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DQB1

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

#hlaDQB1
hlaDQB1<-data.frame(matrix(unlist(hlasplit[7]), byrow=T, nrow=length(hlasplit[[7]])), stringsAsFactors = FALSE)
colnames(hlaDQB1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

template<-hlaDQB1
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,13)))

##UTRS
DQB1_UTR_template<-template
DQB1_UTR_template[,unlist(atlas[7,2])] <- hlaDQB1[unlist(atlas[7,2])]

#CORE EXONS
DQB1_CE_template<-template
DQB1_CE_template[,unlist(atlas[7,3])] <- hlaDQB1[unlist(atlas[7,3])]

#CYTO EXONS
DQB1_CYTO_template<-template
DQB1_CYTO_template[,unlist(atlas[7,4])] <- hlaDQB1[unlist(atlas[7,4])]

#LEADER PEPTIDE
DQB1_LPEP_template<-template
DQB1_LPEP_template[,unlist(atlas[7,5])] <- hlaDQB1[unlist(atlas[7,5])]

#C-DOMAIN
DQB1_CDOM_template<-template
DQB1_CDOM_template[,unlist(atlas[7,6])] <- hlaDQB1[unlist(atlas[7,6])]

#INTRONS
DQB1_INTRONS_template<-template
DQB1_INTRONS_template[,unlist(atlas[7,7])] <- hlaDQB1[unlist(atlas[7,7])]

#TMEXONS
DQB1_TMEXONS_template<-template
DQB1_TMEXONS_template[,unlist(atlas[7,8])] <- hlaDQB1[unlist(atlas[7,8])]


#NON-CORE
DQB1_NONCORE_template<-template
DQB1_NONCORE_template[,unlist(atlas[7,9])] <- hlaDQB1[unlist(atlas[7,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}


