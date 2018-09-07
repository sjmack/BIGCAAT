
setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DRB4

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))


hlaDRB4<-data.frame(matrix(unlist(hlasplit[10]), byrow=T, nrow=length(hlasplit[[10]])), stringsAsFactors = FALSE)
colnames(hlaDRB4)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")


template<-hlaDRB4
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,13)))

##UTRS
DRB4_UTR_template<-template
DRB4_UTR_template[,unlist(atlas[10,2])] <- hlaDRB4[unlist(atlas[10,2])]

#CORE EXONS
DRB4_CE_template<-template
DRB4_CE_template[,unlist(atlas[10,3])] <- hlaDRB4[unlist(atlas[10,3])]

#CYTO EXONS
DRB4_CYTO_template<-template
DRB4_CYTO_template[,unlist(atlas[10,4])] <- hlaDRB4[unlist(atlas[10,4])]

#LEADER PEPTIDE
DRB4_LPEP_template<-template
DRB4_LPEP_template[,unlist(atlas[10,5])] <- hlaDRB4[unlist(atlas[10,5])]

#C-DOMAIN
DRB4_CDOM_template<-template
DRB4_CDOM_template[,unlist(atlas[10,6])] <- hlaDRB4[unlist(atlas[10,6])]

#INTRONS
DRB4_INTRONS_template<-template
DRB4_INTRONS_template[,unlist(atlas[10,7])] <- hlaDRB4[unlist(atlas[10,7])]

#TMEXONS
DRB4_TMEXONS_template<-template
DRB4_TMEXONS_template[,unlist(atlas[10,8])] <- hlaDRB4[unlist(atlas[10,8])]


#NON-CORE
DRB4_NONCORE_template<-template
DRB4_NONCORE_template[,unlist(atlas[10,9])] <- hlaDRB4[unlist(atlas[10,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}


