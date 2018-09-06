setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-B

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

template<-hlaB
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,15)))

##UTRS
B_UTR_template<-template
B_UTR_template[,unlist(atlas[2,2])] <- hlaB[unlist(atlas[2,2])]

#CORE EXONS
B_CE_template<-template
B_CE_template[,unlist(atlas[2,3])] <- hlaB[unlist(atlas[2,3])]

#CYTO EXONS
B_CYTO_template<-template
B_CYTO_template[,unlist(atlas[2,4])] <- hlaB[unlist(atlas[2,4])]

#LEADER PEPTIDE
B_LPEP_template<-template
B_LPEP_template[,unlist(atlas[2,5])] <- hlaB[unlist(atlas[2,5])]

#C-DOMAIN
B_CDOM_template<-template
B_CDOM_template[,unlist(atlas[2,6])] <- hlaB[unlist(atlas[2,6])]

#INTRONS
B_INTRONS_template<-template
B_INTRONS_template[,unlist(atlas[2,7])] <- hlaB[unlist(atlas[2,7])]

#TMEXONS
B_TMEXONS_template<-template
B_TMEXONS_template[,unlist(atlas[2,8])] <- hlaB[unlist(atlas[2,8])]


#NON-CORE
B_NONCORE_template<-template
B_NONCORE_template[,unlist(atlas[2,7])] <- hlaA[unlist(atlas[2,7])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}