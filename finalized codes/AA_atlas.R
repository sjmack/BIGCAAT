#AA determination for exons
#atlas for HLA-A, B, C, DPA1, DPB1, DQB1, and DRB1,
#v 0.8 - 01/15/19
#By:Livia Tran 

###The following script is a set of exon boundaries for HLA-A, B, C, DRB1, DQB1, and DPB1
#These boundaries were calculated with the 3.34.0 IPD-IMGT/HLA alignment files from
#https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments

AA_atlas<-list()
loci<-c("A", "B", "C", "DRB1", "DQB1", "DPB1")

#A
i=1:7
AA_atlas[[loci[1]]] <- data.frame(exon=c(paste("exon", paste(i, i+1, sep=":"), sep="_")), stringsAsFactors = FALSE)
AA_atlas[[loci[[1]]]]$bound <-list(25, 115, 207, 299, 338, 349, 365)

#B
i=1:6
AA_atlas[[loci[2]]] <- data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[2]]]]$bound<-list(25, 115, 207, 299, 339, 350)

#C
i=1:7
AA_atlas[[loci[[3]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[3]]]]$bound<-list(25, 115, 206, 297, 313, 348, 364)

#DRB1
i=1:5
AA_atlas[[loci[[4]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[4]]]]$bound<-list(34, 125, 218, 255, 263)

#DQB1
i=1:5
AA_atlas[[loci[[5]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[5]]]]$bound<-list(34, 124, 218, 255, 267)
  
#DPB1                                
i=1:4
AA_atlas[[loci[[6]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[6]]]]$bound<-list(34, 122, 216, 253)


save(AA_atlas, file="AA_atlas.rda")
                                  
