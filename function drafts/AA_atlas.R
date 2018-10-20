#AA determination for exons
#atlas for HLA-A and HLA-B
#v 0.3 10/16/18
#By:Livia Tran 

AA_atlas<-list()
loci<-c("A", "B")

AA_atlas[[loci[1]]] <- data.frame(exon=paste("exon", seq(1:8)), stringsAsFactors = FALSE)
AA_atlas[[loci[[1]]]]$start <-list(1, 25,115, 207, 299, 338, 349, 365)
AA_atlas[[loci[[1]]]]$end<- list(24, 114, 206, 298, 337, 348, 364, 365)

AA_atlas[[loci[2]]] <- data.frame(exon=paste("exon", seq(1:7)), stringsAsFactors = FALSE)
AA_atlas[[loci[[2]]]]$start<-list(1, 24, 114, 206, 299, 339, 249)
AA_atlas[[loci[[2]]]]$end<-list(23, 113, 205, 298, 337, 348, 363)


save (AA_atlas, file="AA_atlas.rda")


