
setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DRB3

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))


hlaDRB3<-data.frame(matrix(unlist(hlasplit[9]), byrow=T, nrow=length(hlasplit[[9]])), stringsAsFactors = FALSE)
colnames(hlaDRB3)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

template<-hlaDRB3
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,13)))

##UTRS
DRB3_UTR_template<-template
DRB3_UTR_template[,unlist(atlas[9,2])] <- hlaDRB3[unlist(atlas[9,2])]

#CORE EXONS
DRB3_CE_template<-template
DRB3_CE_template[,unlist(atlas[9,3])] <- hlaDRB3[unlist(atlas[9,3])]

#CYTO EXONS
DRB3_CYTO_template<-template
DRB3_CYTO_template[,unlist(atlas[9,4])] <- hlaDRB3[unlist(atlas[9,4])]

#LEADER PEPTIDE
DRB3_LPEP_template<-template
DRB3_LPEP_template[,unlist(atlas[9,5])] <- hlaDRB3[unlist(atlas[9,5])]

#C-DOMAIN
DRB3_CDOM_template<-template
DRB3_CDOM_template[,unlist(atlas[9,6])] <- hlaDRB3[unlist(atlas[9,6])]

#INTRONS
DRB3_INTRONS_template<-template
DRB3_INTRONS_template[,unlist(atlas[9,7])] <- hlaDRB3[unlist(atlas[9,7])]

#TMEXONS
DRB3_TMEXONS_template<-template
DRB3_TMEXONS_template[,unlist(atlas[9,8])] <- hlaDRB3[unlist(atlas[9,8])]


#NON-CORE
DRB3_NONCORE_template<-template
DRB3_NONCORE_template[,unlist(atlas[9,9])] <- hlaDRB3[unlist(atlas[9,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}


