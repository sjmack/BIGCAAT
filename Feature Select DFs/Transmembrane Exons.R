####Feature Sequence Isolation for Transmembrane Exons GFEs


require(stringr)

##see atlas.r, which contains positions for hand-picked sequence features
load("atlas.rda")

##function for reading multiple HLA files in, but keeping them as separate dataframes
multiFileread <- function(filepath, columnnames, skip, clip){
  filenames=list.files(path=filepath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=FALSE, stringsAsFactors=FALSE, skip=skip, col.names = columnnames)})
  datalist<- lapply(datalist, function(x) head(x,-clip)
  )  
  return(datalist)}

hladf<-multiFileread("/Users/liviatran/Desktop/ltmasterscoding/HLA", c("allelename", "gfe"), skip=3, clip=1)
hlasplit<-lapply(hladf, function(x) strsplit(t(as.data.frame(matrix(unlist(strsplit(x[,2],"w")), nrow=nrow(x), byrow=T), stringsAsFactors=FALSE)) [2,],"-"))

#hlaA
hlaAdf<-data.frame(matrix(unlist(hlasplit[1]), byrow=T, nrow=length(hlasplit[[1]])), stringsAsFactors = FALSE)
colnames(hlaAdf)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                    "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR") 

#hlaB
hlaBdf<-data.frame(matrix(unlist(hlasplit[2]), byrow=T, nrow=length(hlasplit[[2]])), stringsAsFactors = FALSE)
colnames(hlaBdf)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                    "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7","3'UTR")


#hlaC
hlaCdf<-data.frame(matrix(unlist(hlasplit[3]), byrow=T, nrow=length(hlasplit[[3]])), stringsAsFactors = FALSE)
colnames(hlaCdf)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                    "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDPA1
hlaDPA1df<-data.frame(matrix(unlist(hlasplit[4]), byrow=T, nrow=length(hlasplit[[4]])), stringsAsFactors = FALSE)
colnames(hlaDPA1df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")

#hlaDPB1
hlaDPB1df<-data.frame(matrix(unlist(hlasplit[5]), byrow=T, nrow=length(hlasplit[[5]])), stringsAsFactors = FALSE)
colnames(hlaDPB1df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                       "Intron 4", "Exon 5","3'UTR")
#hlaDQA1
hlaDQA1df<-data.frame(matrix(unlist(hlasplit[6]), byrow=T, nrow=length(hlasplit[[6]])), stringsAsFactors = FALSE)
colnames(hlaDQA1df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")

#hlaDQB1
hlaDQB1df<-data.frame(matrix(unlist(hlasplit[7]), byrow=T, nrow=length(hlasplit[[7]])), stringsAsFactors = FALSE)
colnames(hlaDQB1df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                       "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB1
hlaDRB1df<-data.frame(matrix(unlist(hlasplit[8]), byrow=T, nrow=length(hlasplit[[8]])), stringsAsFactors = FALSE)
colnames(hlaDRB1df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                       "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB3
hlaDRB3df<-data.frame(matrix(unlist(hlasplit[9]), byrow=T, nrow=length(hlasplit[[9]])), stringsAsFactors = FALSE)
colnames(hlaDRB3df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                       "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB4
hlaDRB4df<-data.frame(matrix(unlist(hlasplit[10]), byrow=T, nrow=length(hlasplit[[10]])), stringsAsFactors = FALSE)
colnames(hlaDRB4df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                       "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")

#hlaDRB5
hlaDRB5df<-data.frame(matrix(unlist(hlasplit[11]), byrow=T, nrow=length(hlasplit[[11]])), stringsAsFactors = FALSE)
colnames(hlaDRB5df)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                       "Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")  

##prepares for building a list for UTR selected features for each loci
loci<-c("HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5")
TMexons<-list()

#function to make zeroed out templates based on each loci - then fills in with appropriate information
#based on coordinates of the atlas 
featureselect<-function(dataframe, x, y, GFEsplitDF){
  dataframe[,1:ncol(dataframe)] <- as.data.frame(rbind(rep(0,length(dataframe))))
  dataframe[,unlist(atlas[x,y])] <- GFEsplitDF[unlist(atlas[x,y])]
  return(dataframe)}

TMexons[[loci[1]]]<-featureselect(hlaAdf, 1,8, hlaAdf)
TMexons[[loci[2]]]<-featureselect(hlaBdf, 2,8, hlaBdf)
TMexons[[loci[3]]]<-featureselect(hlaCdf, 3,8, hlaCdf)
TMexons[[loci[4]]]<-featureselect(hlaDPA1df, 4,8, hlaDPA1df)
TMexons[[loci[5]]]<-featureselect(hlaDPB1df, 5,8, hlaDPB1df)
TMexons[[loci[6]]]<-featureselect(hlaDQA1df, 6,8, hlaDQA1df)
TMexons[[loci[7]]]<-featureselect(hlaDQB1df, 7,8, hlaDQB1df)
TMexons[[loci[8]]]<-featureselect(hlaDRB1df, 8,8, hlaDRB1df)
TMexons[[loci[9]]]<-featureselect(hlaDRB3df, 9,8, hlaDRB3df)
TMexons[[loci[10]]]<-featureselect(hlaDRB4df, 10, 8, hlaDRB4df)
TMexons[[loci[11]]]<-featureselect(hlaDRB5df, 11,8, hlaDRB5df)

GFEpaster<-function(seqfeaturelist, seqfeatureposition, lociDF, BSGdf, BSGdfposition, zerolocus){
  seqfeaturelist[[seqfeatureposition]][nrow(seqfeaturelist[[seqfeatureposition]]) + 1,] = max(as.numeric(as.character(unlist(lociDF[,2:ncol(lociDF)]))))+1
  cols<- 1:length(seqfeaturelist[[seqfeatureposition]])
  GFEs<-str_replace_all(paste(gsub("[*].*", "", unlist(BSGdf[[BSGdfposition]][1])),paste("w", apply(seqfeaturelist[[seqfeatureposition]][,cols] , 1 , paste, collapse="-"))), " ", "")
  BSGdf[[BSGdfposition]][nrow(BSGdf[[BSGdfposition]]) + 1,] = list(zerolocus, NA)
  featureselectDF<-cbind(BSGdf[[BSGdfposition]][[1]], GFEs)
  return(featureselectDF)}



#builds a new list based on previously defined loci -- list contains allele names and sequence feature specific GFE notations
TMexonsIsolated<-list()
TMexonsIsolated[[loci[1]]]<-GFEpaster(TMexons, 1, hlaAdf, hladf, 1, "HLA-A*00:00")
TMexonsIsolated[[loci[2]]]<-GFEpaster(TMexons, 2, hlaBdf, hladf, 2, "HLA-B*00:00")
TMexonsIsolated[[loci[3]]]<-GFEpaster(TMexons, 3, hlaCdf, hladf, 3, "HLA-C*00:00")
TMexonsIsolated[[loci[4]]]<-GFEpaster(TMexons, 4, hlaDPA1df, hladf, 4, "HLA-DPA1*00:00")
TMexonsIsolated[[loci[5]]]<-GFEpaster(TMexons, 5, hlaDPB1df, hladf, 5, "HLA-DPB1*00:00")
TMexonsIsolated[[loci[6]]]<-GFEpaster(TMexons, 6, hlaDQA1df, hladf, 6, "HLA-DQA1*00:00")
TMexonsIsolated[[loci[7]]]<-GFEpaster(TMexons, 7, hlaDQB1df, hladf, 7, "HLA-DQB1*00:00")
TMexonsIsolated[[loci[8]]]<-GFEpaster(TMexons, 8, hlaDRB1df, hladf, 8, "HLA-DRB1*00:00")
TMexonsIsolated[[loci[9]]]<-GFEpaster(TMexons, 9, hlaDRB3df, hladf, 9, "HLA-DRB3*00:00")
TMexonsIsolated[[loci[10]]]<-GFEpaster(TMexons, 10, hlaDRB4df, hladf, 10, "HLA-DRB4*00:00")
TMexonsIsolated[[loci[11]]]<-GFEpaster(TMexons, 11, hlaDRB5df, hladf, 11, "HLA-DRB5*00:00")  

