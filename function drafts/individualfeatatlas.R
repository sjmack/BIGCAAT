individualfeatatlas<-data.frame(locus=names(framework))

individualfeatatlas$exon1<-list(2, 2, 2, 2,
                                2, 2, 2, 2,
                                2, 2, 2)
                                

individualfeatatlas$intron1<-list(3,3,3,3,
                                  3,3,3,3,
                                  3,3,3)

individualfeatatlas$exon2<-list(4,4,4,4,
                                4,4,4,4,
                                4,4,4)
                                
individualfeatatlas$intron2<-list(5, 5, 5, 5,
                                  5, 5, 5, 5,
                                  5, 5,5)

individualfeatatlas$exon3<-list(6, 6, 6, 6, 
                                6, 6, 6, 6,
                                6, 6, 6)

individualfeatatlas$intron3<-list( 7, 7, 7, 7,
                                   7, 7, 7, 7,
                                   7, 7, 7)


individualfeatatlas$exon4<-list(8, 8, 8, 8,
                                8, 8, 8, 8,
                                8,8,8)

individualfeatatlas$intron4<-list(9, 9, 9, 0,
                                  9, 0, 9, 9,
                                  9,9,9)

individualfeatatlas$exon5<-list(10, 10, 10, 0,
                                10, 0, 10, 10,
                                10,10,10)

individualfeatatlas$intron5<-list(11, 11, 11, 0,
                                  0, 0, 11, 11,
                                  11, 11, 11)

individualfeatatlas$exon6<-list(12, 12, 12, 0,
                                  0, 0, 12, 12,
                                  12, 12, 12)

individualfeatatlas$intron6<-list(13, 13, 13, 0,
                                0, 0, 0, 0,
                                0, 0, 0)

individualfeatatlas$exon7<-list(14, 14, 14, 0,
                                  0, 0, 0, 0,
                                  0, 0, 0)
individualfeatatlas$intron7<-list(15, 0, 15, 0,
                                0, 0, 0, 0,
                                0, 0, 0)

individualfeatatlas$exon8<-list(16, 0, 16, 0,
                                  0, 0, 0, 0,
                                  0, 0, 0)
setwd("~/Desktop/ltmasterscoding")
save(individualfeatatlas, file="individualfeatatlas.rda")
