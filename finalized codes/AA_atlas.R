#AA determination for exons
#atlas for HLA-A, B, C, DPA1, DPB1, DQB1, and DRB1,
#v 0.7 - 12/18/18
#By:Livia Tran 

###The following script is a set of exon boundaries for HLA-A, B, C, DRB1, DQB1, and DPB1
#These boundaries were calculated with the 3.34.0 IPD-IMGT/HLA alignment files from
#https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments

AA_atlas<-list()
loci<-c("A", "B", "C", "DRB1", "DQB1", "DPB1")

#A
i=1:7
AA_atlas[[loci[1]]] <- data.frame(exon=c(paste("exon", paste(i, i+1, sep=":"), sep="_")), stringsAsFactors = FALSE)
AA_atlas[[loci[[1]]]]$bound <-list(1, 91, 183, 275, 314, 325, 341)

#B
i=1:6
AA_atlas[[loci[2]]] <- data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[2]]]]$bound<-list(1, 91, 183, 275, 315, 326)

#C
i=1:7
AA_atlas[[loci[[3]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[3]]]]$bound<-list(1, 91, 182, 273, 313, 324, 340)

#DRB1
i=1:5
AA_atlas[[loci[[4]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[4]]]]$bound<-list(5, 96, 189, 226, 234)

#DQB1
i=1:5
AA_atlas[[loci[[5]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[5]]]]$bound<-list(5, 95, 189, 226, 235)
  
#DPB1                                
i=1:4
AA_atlas[[loci[[6]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[6]]]]$bound<-list(5, 93, 187, 224)


save(AA_atlas, file="AA_atlas.rda")
                                  
