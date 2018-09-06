setwd("~/Desktop/ltmasterscoding")

#################################Feature Selection for HLA-A

###see atlas.R

########FOR HLA-A
hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/BSG 3.4", c("allelename", "gfe"), skip=3, clip=1)

#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))


#turns split list of individual gene features into matrix, and then to dataframe with every
#gene feature getting its own column -- names features 

#hlaA
hlaA<-data.frame(matrix(unlist(hlasplit[1]), byrow=T, nrow=length(hlasplit[[1]])), stringsAsFactors = FALSE)
colnames(hlaA)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                  "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR") 


template<-hlaA
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,17)))

##UTRS
A_UTR_template<-template
A_UTR_template[,unlist(atlas[1,2])] <- hlaA[unlist(atlas[1,2])]

#CORE EXONS
A_CE_template<-template
A_CE_template[,unlist(atlas[1,3])] <- hlaA[unlist(atlas[1,3])]

#CYTO EXONS
A_CYTO_template<-template
A_CYTO_template[,unlist(atlas[1,4])] <- hlaA[unlist(atlas[1,4])]

#LEADER PEPTIDE
A_LPEP_template<-template
A_LPEP_template[,unlist(atlas[1,5])] <- hlaA[unlist(atlas[1,5])]

#C-DOMAIN
A_CDOM_template<-template
A_CDOM_template[,unlist(atlas[1,6])] <- hlaA[unlist(atlas[1,6])]

#INTRONS
A_INTRONS_template<-template
A_INTRONS_template[,unlist(atlas[1,7])] <- hlaA[unlist(atlas[1,7])]

#TMEXONS
A_TMEXONS_template<-template
A_TMEXONS_template[,unlist(atlas[1,8])] <- hlaA[unlist(atlas[1,8])]


#NON-CORE
A_NONCORE_template<-template
A_NONCORE_template[,unlist(atlas[1,7])] <- hlaA[unlist(atlas[1,7])]


#for reading multiple files into one dataframe
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}