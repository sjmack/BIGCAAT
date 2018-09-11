####Feature Sequence Isolation for Intron GFEs


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
IntronList<-list()

#function to make zeroed out templates based on each loci - then fills in with appropriate information
#based on coordinates of the atlas 
featureselect<-function(dataframe, x, y, GFEsplitDF){
  dataframe[,1:ncol(dataframe)] <- as.data.frame(rbind(rep(0,length(dataframe))))
  dataframe[,unlist(atlas[x,y])] <- GFEsplitDF[unlist(atlas[x,y])]
  return(dataframe)
  
  IntronList[[loci[1]]]<-featureselect(hlaAdf, 1,7, hlaAdf)
  IntronList[[loci[2]]]<-featureselect(hlaBdf, 2,7, hlaBdf)
  IntronList[[loci[3]]]<-featureselect(hlaCdf, 3,7, hlaCdf)
  IntronList[[loci[4]]]<-featureselect(hlaDPA1df, 4,7, hlaDPA1df)
  IntronList[[loci[5]]]<-featureselect(hlaDPB1df, 5,7, hlaDPB1df)
  IntronList[[loci[6]]]<-featureselect(hlaDQA1df, 6,7, hlaDQA1df)
  IntronList[[loci[7]]]<-featureselect(hlaDQB1df, 7,7, hlaDQB1df)
  IntronList[[loci[8]]]<-featureselect(hlaDRB1df, 8,7, hlaDRB1df)
  IntronList[[loci[9]]]<-featureselect(hlaDRB3df, 9,7, hlaDRB3df)
  IntronList[[loci[10]]]<-featureselect(hlaDRB4df, 10, 7, hlaDRB4df)
  IntronList[[loci[11]]]<-featureselect(hlaDRB5df, 11,7, hlaDRB5df)
  
  GFEpaster<-function(seqfeaturelist, seqfeatureposition, lociDF, BSGdf, BSGdfposition, zerolocus){
    seqfeaturelist[[seqfeatureposition]][nrow(seqfeaturelist[[seqfeatureposition]]) + 1,] = max(as.numeric(as.character(unlist(lociDF[,2:ncol(lociDF)]))))+1
    cols<- 1:length(seqfeaturelist[[seqfeatureposition]])
    GFEs<-str_replace_all(paste(gsub("[*].*", "", unlist(BSGdf[[BSGdfposition]][1])),paste("w", apply(seqfeaturelist[[seqfeatureposition]][,cols] , 1 , paste, collapse="-"))), " ", "")
    BSGdf[[BSGdfposition]][nrow(BSGdf[[BSGdfposition]]) + 1,] = list(zerolocus, NA)
    featureselectDF<-cbind(BSGdf[[BSGdfposition]][[1]], GFEs)
    return(featureselectDF)}
  
  
  
  #builds a new list based on previously defined loci -- list contains allele names and sequence feature specific GFE notations
  IntronIsolated<-list()
  IntronIsolated[[loci[1]]]<-GFEpaster(IntronList, 1, hlaAdf, hladf, 1, "HLA-A*00:00")
  IntronIsolated[[loci[2]]]<-GFEpaster(IntronList, 2, hlaBdf, hladf, 2, "HLA-B*00:00")
  IntronIsolated[[loci[3]]]<-GFEpaster(IntronList, 3, hlaCdf, hladf, 3, "HLA-C*00:00")
  IntronIsolated[[loci[4]]]<-GFEpaster(IntronList, 4, hlaDPA1df, hladf, 4, "HLA-DPA1*00:00")
  IntronIsolated[[loci[5]]]<-GFEpaster(IntronList, 5, hlaDPB1df, hladf, 5, "HLA-DPB1*00:00")
  IntronIsolated[[loci[6]]]<-GFEpaster(IntronList, 6, hlaDQA1df, hladf, 6, "HLA-DQA1*00:00")
  IntronIsolated[[loci[7]]]<-GFEpaster(IntronList, 7, hlaDQB1df, hladf, 7, "HLA-DQB1*00:00")
  IntronIsolated[[loci[8]]]<-GFEpaster(IntronList, 8, hlaDRB1df, hladf, 8, "HLA-DRB1*00:00")
  IntronIsolated[[loci[9]]]<-GFEpaster(IntronList, 9, hlaDRB3df, hladf, 9, "HLA-DRB3*00:00")
  IntronIsolated[[loci[10]]]<-GFEpaster(IntronList, 10, hlaDRB4df, hladf, 10, "HLA-DRB4*00:00")
  IntronIsolated[[loci[11]]]<-GFEpaster(IntronList, 11, hlaDRB5df, hladf, 11, "HLA-DRB5*00:00")  
  
  