
setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DPB1

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

#hlaDQA1
hlaDQA1<-data.frame(matrix(unlist(hlasplit[6]), byrow=T, nrow=length(hlasplit[[6]])), stringsAsFactors = FALSE)
colnames(hlaDQA1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")

template<-hlaDQA1
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,9)))

##UTRS
DQA1_UTR_template<-template
DQA1_UTR_template[,unlist(atlas[6,2])] <- hlaDQA1[unlist(atlas[6,2])]

#CORE EXONS
DQA1_CE_template<-template
DQA1_CE_template[,unlist(atlas[6,3])] <- hlaDQA1[unlist(atlas[6,3])]

#CYTO EXONS
DQA1_CYTO_template<-template
DQA1_CYTO_template[,unlist(atlas[6,4])] <- hlaDQA1[unlist(atlas[6,4])]

#LEADER PEPTIDE
DQA1_LPEP_template<-template
DQA1_LPEP_template[,unlist(atlas[6,5])] <- hlaDQA1[unlist(atlas[6,5])]

#C-DOMAIN
DQA1_CDOM_template<-template
DQA1_CDOM_template[,unlist(atlas[6,6])] <- hlaDQA1[unlist(atlas[6,6])]

#INTRONS
DQA1_INTRONS_template<-template
DQA1_INTRONS_template[,unlist(atlas[6,7])] <- hlaDQA1[unlist(atlas[6,7])]

#TMEXONS
DQA1_TMEXONS_template<-template
DQA1_TMEXONS_template[,unlist(atlas[6,8])] <- hlaDQA1[unlist(atlas[6,8])]


#NON-CORE
DQA1_NONCORE_template<-template
DQA1_NONCORE_template[,unlist(atlas[6,9])] <- hlaDQA1[unlist(atlas[6,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}


