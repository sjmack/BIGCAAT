#AA determination for exons
#atlas for HLA-A
#v 0.1
#By:Livia Tran 


#AA atlas formation 
AA_atlas <- data.frame(exon=paste("exon", seq(1:8)), stringsAsFactors = FALSE)

AA_atlas$start <- paste("position", list(1, 25,115, 207, 299, 338, 349, 365), sep=".")

AA_atlas$end<- paste("position", list(24, 114, 206, 298, 337, 348, 364, 365), sep=".")

