library(BIGDAWG)

View(HLA_data)
hla_data<-BDgenotypeconversion("/Users/liviatran/Desktop/ltmasterscoding/HLA", "samplebdHLA.csv", "Allelelist.3310.txt", CWDverify())

omittedata<-list()

for(i in 3:ncol(hla_data)){
omittedata[[i]]<-na.omit(hla_data[,i])}



hla_data<-BDgenotypeconversion("/Users/liviatran/Desktop/ltmasterscoding/HLA", "samplebdHLA.csv", "Allelelist.3310.txt", CWDverify())
hla_data[is.na(hla_data)] <- ""
View(na.omit(hla_data))
hlasplit<-list()

for(i in 3:ncol(hla_data)){
hlasplit[[i]]<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hla_data[,i],"w")), nrow=nrow(hla_data), byrow=T), stringsAsFactors=FALSE)) [2,], "-")}


hlaAdf<-data.frame(matrix(unlist(hlasplit[[3]]), byrow=T, nrow=length(hlasplit[[3]])), stringsAsFactors = FALSE)

colnames(hlaAdf)<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4",
                  "Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR") 


template<-hlaAdf
template[,1:ncol(template)] <- as.data.frame(rbind(rep(0,17)))
A_UTR_template<-template
A_UTR_template[,unlist(atlas[1,2])] <- hlaAdf[unlist(atlas[1,2])]

newdf<- data.frame("GFE"=str_replace_all(paste("HLA-Aw", apply( A_UTR_template[ , 1:17] , 1 , paste , collapse = "-" )), " ", ""))
for(i in 3:ncol(newdf)){
  newdf[i]<-paste(sapply(hlasplit[[i]], "[", 1), sapply(hlasplit[[i]], "[", 2), sep=":")
  newdf[,i]<-paste(colnames(newdf[i]),newdf[,i],sep="*")}
