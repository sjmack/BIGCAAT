
setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-DPB1

###see atlas.R

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

hlaDPB1<-data.frame(matrix(unlist(hlasplit[5]), byrow=T, nrow=length(hlasplit[[5]])), stringsAsFactors = FALSE)
colnames(hlaDPB1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5","3'UTR")
template<-hlaDPB1
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,11)))

##UTRS
DPB1_UTR_template<-template
DPB1_UTR_template[,unlist(atlas[5,2])] <- hlaDPB1[unlist(atlas[5,2])]

#CORE EXONS
DPB1_CE_template<-template
DPB1_CE_template[,unlist(atlas[5,3])] <- hlaDPB1[unlist(atlas[5,3])]

#CYTO EXONS
DPB1_CYTO_template<-template
DPB1_CYTO_template[,unlist(atlas[5,4])] <- hlaDPB1[unlist(atlas[5,4])]

#LEADER PEPTIDE
DPB1_LPEP_template<-template
DPB1_LPEP_template[,unlist(atlas[5,5])] <- hlaDPB1[unlist(atlas[5,5])]

#C-DOMAIN
DPB1_CDOM_template<-template
DPB1_CDOM_template[,unlist(atlas[5,6])] <- hlaDPB1[unlist(atlas[5,6])]

#INTRONS
DPB1_INTRONS_template<-template
DPB1_INTRONS_template[,unlist(atlas[5,7])] <- hlaDPB1[unlist(atlas[5,7])]

#TMEXONS
DPB1_TMEXONS_template<-template
DPB1_TMEXONS_template[,unlist(atlas[5,8])] <- hlaDPB1[unlist(atlas[5,8])]


#NON-CORE
DPB1_NONCORE_template<-template
DPB1_NONCORE_template[,unlist(atlas[5,9])] <- hlaDPB1[unlist(atlas[5,9])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}


