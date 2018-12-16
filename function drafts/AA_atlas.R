#AA determination for exons
#atlas for HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1, DRB3/4/5
#v 0.6 - 10/27/18
#By:Livia Tran 

AA_atlas<-list()
loci<-c("A", "B", "C", "DRB", "DPB1", "DQB1") 

i=1:7
AA_atlas[[loci[1]]] <- data.frame(exon=c(paste("exon", paste(i, i+1, sep=":"), sep="_")), stringsAsFactors = FALSE)
AA_atlas[[loci[[1]]]]$bound <-list(1, 91, 183, 275, 314, 325, 341)

i=1:6
AA_atlas[[loci[2]]] <- data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[2]]]]$bound<-list(1, 91, 183, 275, 315, 326)


i=1:7
AA_atlas[[loci[[3]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[3]]]]$bound<-list(1, 91, 182, 273, 313, 324, 340)

i=1:5
AA_atlas[[loci[[4]]]]<-data.frame(exon=paste("exon", paste(i, i+1, sep=":"), sep="_"), stringsAsFactors = FALSE)
AA_atlas[[loci[[4]]]]$bound<-list(5, 95, 189, 226, 234)

View(AA_atlas)
