
setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DRB1

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))


hlaDRB1<-data.frame(matrix(unlist(hlasplit[8]), byrow=T, nrow=length(hlasplit[[8]])), stringsAsFactors = FALSE)
colnames(hlaDRB1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

template<-hlaDRB1
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,13)))

##UTRS
DRB1_UTR_template<-template
DRB1_UTR_template[,unlist(atlas[8,2])] <- hlaDRB1[unlist(atlas[8,2])]

#CORE EXONS
DRB1_CE_template<-template
DRB1_CE_template[,unlist(atlas[8,3])] <- hlaDRB1[unlist(atlas[8,3])]

#CYTO EXONS
DRB1_CYTO_template<-template
DRB1_CYTO_template[,unlist(atlas[8,4])] <- hlaDRB1[unlist(atlas[8,4])]

#LEADER PEPTIDE
DRB1_LPEP_template<-template
DRB1_LPEP_template[,unlist(atlas[8,5])] <- hlaDRB1[unlist(atlas[8,5])]

#C-DOMAIN
DRB1_CDOM_template<-template
DRB1_CDOM_template[,unlist(atlas[8,6])] <- hlaDRB1[unlist(atlas[8,6])]

#INTRONS
DRB1_INTRONS_template<-template
DRB1_INTRONS_template[,unlist(atlas[8,7])] <- hlaDRB1[unlist(atlas[8,7])]

#TMEXONS
DRB1_TMEXONS_template<-template
DRB1_TMEXONS_template[,unlist(atlas[8,8])] <- hlaDRB1[unlist(atlas[8,8])]


#NON-CORE
DRB1_NONCORE_template<-template
DRB1_NONCORE_template[,unlist(atlas[8,9])] <- hlaDRB1[unlist(atlas[8,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}


