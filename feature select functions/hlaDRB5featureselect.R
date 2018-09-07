
setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DRB5

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))


#hlaDRB5
hlaDRB5<-data.frame(matrix(unlist(hlasplit[11]), byrow=T, nrow=length(hlasplit[[11]])), stringsAsFactors = FALSE)
colnames(hlaDRB5)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")  

template<-hlaDRB5
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,13)))

##UTRS
DRB5_UTR_template<-template
DRB5_UTR_template[,unlist(atlas[11,2])] <- hlaDRB5[unlist(atlas[11,2])]

#CORE EXONS
DRB5_CE_template<-template
DRB5_CE_template[,unlist(atlas[11,3])] <- hlaDRB5[unlist(atlas[11,3])]

#CYTO EXONS
DRB5_CYTO_template<-template
DRB5_CYTO_template[,unlist(atlas[11,4])] <- hlaDRB5[unlist(atlas[11,4])]

#LEADER PEPTIDE
DRB5_LPEP_template<-template
DRB5_LPEP_template[,unlist(atlas[11,5])] <- hlaDRB5[unlist(atlas[11,5])]

#C-DOMAIN
DRB5_CDOM_template<-template
DRB5_CDOM_template[,unlist(atlas[11,6])] <- hlaDRB5[unlist(atlas[11,6])]

#INTRONS
DRB5_INTRONS_template<-template
DRB5_INTRONS_template[,unlist(atlas[11,7])] <- hlaDRB5[unlist(atlas[11,7])]

#TMEXONS
DRB5_TMEXONS_template<-template
DRB5_TMEXONS_template[,unlist(atlas[11,8])] <- hlaDRB5[unlist(atlas[11,8])]


#NON-CORE
DRB5_NONCORE_template<-template
DRB5_NONCORE_template[,unlist(atlas[11,9])] <- hlaDRB5[unlist(atlas[11,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}


