setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DPA1

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

hlaDPA1<-data.frame(matrix(unlist(hlasplit[4]), byrow=T, nrow=length(hlasplit[[4]])), stringsAsFactors = FALSE)
colnames(hlaDPA1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")


template<-hlaDPA1
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,9)))

##UTRS
DPA1_UTR_template<-template
DPA1_UTR_template[,unlist(atlas[4,2])] <- hlaDPA1[unlist(atlas[4,2])]

#CORE EXONS
DPA1_CE_template<-template
DPA1_CE_template[,unlist(atlas[4,3])] <- hlaDPA1[unlist(atlas[4,3])]

#CYTO EXONS
DPA1_CYTO_template<-template
DPA1_CYTO_template[,unlist(atlas[4,4])] <- hlaDPA1[unlist(atlas[4,4])]

#LEADER PEPTIDE
DPA1_LPEP_template<-template
DPA1_LPEP_template[,unlist(atlas[4,5])] <- hlaDPA1[unlist(atlas[4,5])]

#C-DOMAIN
DPA1_CDOM_template<-template
DPA1_CDOM_template[,unlist(atlas[4,6])] <- hlaDPA1[unlist(atlas[4,6])]

#INTRONS
DPA1_INTRONS_template<-template
DPA1_INTRONS_template[,unlist(atlas[4,7])] <- hlaDPA1[unlist(atlas[4,7])]

#TMEXONS
DPA1_TMEXONS_template<-template
DPA1_TMEXONS_template[,unlist(atlas[4,8])] <- hlaDPA1[unlist(atlas[4,8])]


#NON-CORE
DPA1_NONCORE_template<-template
DPA1_NONCORE_template[,unlist(atlas[4,9])] <- hlaDPA1[unlist(atlas[4,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}