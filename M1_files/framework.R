###framework formation
##version 2.0
#By: Livia Tran & S.J. Mack. 


framework <- list()
locus<-c("HLA-A","HLA-B","HLA-C", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1","HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5")
framework[[locus[1]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8", "3'UTR") 
framework[[locus[2]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7","3'UTR")
framework[[locus[3]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6", "Intron 6", "Exon 7", "Intron 7", "Exon 8","3'UTR")
framework[[locus[4]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")
framework[[locus[5]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5","3'UTR")
framework[[locus[6]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","3'UTR")
framework[[locus[7]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")
framework[[locus[8]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")
framework[[locus[9]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")
framework[[locus[10]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")
framework[[locus[11]]]<-c("5'UTR", "Exon 1", "Intron 1","Exon 2","Intron 2","Exon 3", "Intron 3", "Exon 4","Intron 4", "Exon 5", "Intron 5", "Exon 6","3'UTR")  
save(framework,file="framework.rda")
