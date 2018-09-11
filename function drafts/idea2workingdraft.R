setwd("~/Desktop/ltmasterscoding")

##atlas formation 
atlas <- data.frame(locus=c("HLA-A","HLA-B","HLA-C",
                            "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
                            "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5"),stringsAsFactors = FALSE)
atlas$UTRS <- list(c(1,17),c(1,15),c(1,17),
                   c(1,9), c(1,11), c(1,9), c(1,13), 
                   c(1,13), c(1,13), c(1,13), c(1,13))
atlas$core <- list(c(4,6),c(4,6),c(4,6),
                   4, 4, 4, 4,
                   4, 4, 4, 4)
atlas$cytoexons <- list(c(12,14,16),c(12,14),c(12,14,16),
                        8, 10, 8, c(10, 12), 
                        c(10,12), c(10,12), c(10,12), c(10,12))
atlas$leaderpep<-list(2, 2,2,
                      2,2,2,2,
                      2,2,2,2)
atlas$cdomain<-list(8,8,8,
                    8,8,8,8,
                    8,8,8,8)
atlas$introns<-list(c(3, 5, 7, 9, 11, 13, 15), c(3,5,7,9,11,13),c(3, 5, 7, 9, 11, 13, 15),
                    c(3,5,7), c(3,5,7,9), c(3,5,7), c(3,5,7,9,11),
                    c(3,5,7,9,11),c(3,5,7,9,11),c(3,5,7,9,11),c(3,5,7,9,11))
atlas$tmexons<-list(10,10,10,
                    8,8,8,8,
                    8,8,8,8)
atlas$noncore<-list(c(2, 8, 10, 12, 14, 16), c(2, 8, 10, 12), c(2, 8, 10, 12, 14, 16),
                    c(2, 6, 8), c(2, 6, 8, 10), c(2, 6, 8), c(2, 6, 8, 10, 12),
                    c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12))

#for file merging
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
)  
  return(datalist)}

View(HLA_data)
hladf[[9]]
hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/HLA", c("allelename", "gfe"), skip=3, clip=1)
#splits hladf with each respective gene feature getting its own column
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))


#turns split list of individual gene features into matrix, and then to dataframe with every
#gene feature getting its own column -- names features 

#hlaA
hlaA<-data.frame(matrix(unlist(hlasplit[1]), byrow=T, nrow=length(hlasplit[[1]])), stringsAsFactors = FALSE)
colnames(hlaA)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
            "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR") 


max(as.numeric(as.character(unlist(hlaA[,2:ncol(hlaA)]))))

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
A_NONCORE_template[,unlist(atlas[1,7])] <- hlaA[unlist(atlas[1,9])]

#hlaB
hlaB<-data.frame(matrix(unlist(hlasplit[2]), byrow=T, nrow=length(hlasplit[[2]])), stringsAsFactors = FALSE)
colnames(hlaB)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
  "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7","3'UTR")



#hlaC
hlaC<-data.frame(matrix(unlist(hlasplit[3]), byrow=T, nrow=length(hlasplit[[3]])), stringsAsFactors = FALSE)
colnames(hlaC)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                  "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDPA1
hlaDPA1<-data.frame(matrix(unlist(hlasplit[4]), byrow=T, nrow=length(hlasplit[[4]])), stringsAsFactors = FALSE)
colnames(hlaDPA1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")

#hlaDPB1
hlaDPB1<-data.frame(matrix(unlist(hlasplit[5]), byrow=T, nrow=length(hlasplit[[5]])), stringsAsFactors = FALSE)
colnames(hlaDPB1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                  "Intron 4", "Exon 5","3'UTR")
#hlaDQA1
hlaDQA1<-data.frame(matrix(unlist(hlasplit[6]), byrow=T, nrow=length(hlasplit[[6]])), stringsAsFactors = FALSE)
colnames(hlaDQA1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")

#hlaDQB1
hlaDQB1<-data.frame(matrix(unlist(hlasplit[7]), byrow=T, nrow=length(hlasplit[[7]])), stringsAsFactors = FALSE)
colnames(hlaDQB1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB1
hlaDRB1<-data.frame(matrix(unlist(hlasplit[8]), byrow=T, nrow=length(hlasplit[[8]])), stringsAsFactors = FALSE)
colnames(hlaDRB1)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB3
hlaDRB3<-data.frame(matrix(unlist(hlasplit[9]), byrow=T, nrow=length(hlasplit[[9]])), stringsAsFactors = FALSE)
colnames(hlaDRB3)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB4
hlaDRB4<-data.frame(matrix(unlist(hlasplit[10]), byrow=T, nrow=length(hlasplit[[10]])), stringsAsFactors = FALSE)
colnames(hlaDRB4)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB5
hlaDRB5<-data.frame(matrix(unlist(hlasplit[11]), byrow=T, nrow=length(hlasplit[[11]])), stringsAsFactors = FALSE)
colnames(hlaDRB5)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                     "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")  


max(as.numeric(as.character(unlist(hlaDRB5[,2:ncol(hlaDRB5)]))))
hlaDRB5<-max(as.numeric(as.character(unlist(hlaDRB5[,2:ncol(hlaDRB5)]))))+1 


HLAdata[is.na(HLAdata)] <- "0"

ifelse(HLAdata$DRB5=="NA", )