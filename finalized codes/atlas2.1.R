##atlas
#By: Livia Tran
#September 21, 2019
#v 2.1

atlas<-data.frame(locus=names(framework))

atlas$fiveUTR<-list(1, 1, 1, 1,
                 1, 1, 1, 1,
                 1, 1, 1)


atlas$exon1_lpeptide<-list(2, 2, 2, 2,
                                         2, 2, 2, 2,
                                         2, 2, 2)


atlas$intron1<-list(3,3,3,3,
                                  3,3,3,3,
                                  3,3,3)


atlas$exon2<-list(4,4,4,4,
                                4,4,4,4,
                                4,4,4)


atlas$intron2<-list(5, 5, 5, 5,
                                  5, 5, 5, 5,
                                  5, 5,5)

atlas$exon3<-list(6, 6, 6, 6, 
                                6, 6, 6, 6,
                                6, 6, 6)


atlas$intron3<-list( 7, 7, 7, 7,
                                   7, 7, 7, 7,
                                   7, 7, 7)


atlas$exon4_cdomain<-list(8, 8, 8, 8,
                                8, 8, 8, 8,
                                8,8,8)



atlas$intron4<-list(9, 9, 9, 0,
                                  9, 0, 9, 9,
                                  9,9,9)

atlas$exon5<-list(10, 10, 10, 0,
                                10, 0, 10, 10,
                                10,10,10)

atlas$intron5<-list(11, 11, 11, 0,
                                  0, 0, 11, 11,
                                  11, 11, 11)
atlas$exon6<-list(12, 12, 12, 0,
                                0, 0, 12, 12,
                                12, 12, 12)

atlas$intron6<-list(13, 13, 13, 0,
                                  0, 0, 0, 0,
                                  0, 0, 0)

atlas$exon7<-list(14, 14, 14, 0,
                                0, 0, 0, 0,
                                0, 0, 0)
atlas$intron7<-list(15, 0, 15, 0,
                                  0, 0, 0, 0,
                                  0, 0, 0)

atlas$exon8<-list(16, 0, 16, 0,
                                0, 0, 0, 0,
                                0, 0, 0)

atlas$threeUTR<-list(17, 15, 17, 9, 
                     11, 9, 13, 13,
                     13, 13, 13)


atlas$UTRS <- list(c(1,17),c(1,15),c(1,17),
                   c(1,9), c(1,11), c(1,9), c(1,13), 
                   c(1,13), c(1,13), c(1,13), c(1,13))


atlas$core <- list(c(4,6),c(4,6),c(4,6),
                   4, 4, 4, 4,
                   4, 4, 4, 4)

atlas$cytoexons <- list(c(12,14,16),c(12,14),c(12,14,16),
                        8, 10, 8, c(10, 12), 
                        c(10,12), c(10,12), c(10,12), c(10,12))

atlas$introns<-list(c(3, 5, 7, 9, 11, 13, 15), c(3,5,7,9,11,13),c(3, 5, 7, 9, 11, 13, 15),
                    c(3,5,7), c(3,5,7,9), c(3,5,7), c(3,5,7,9,11),
                    c(3,5,7,9,11),c(3,5,7,9,11),c(3,5,7,9,11),c(3,5,7,9,11))

atlas$tmexons<-list(10,10,10,
                    8,8,8,8,
                    8,8,8,8)
atlas$noncore<-list(c(2, 8, 10, 12, 14, 16), c(2, 8, 10, 12), c(2, 8, 10, 12, 14, 16),
                    c(2, 6, 8), c(2, 6, 8, 10), c(2, 6, 8), c(2, 6, 8, 10, 12),
                    c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12), c(2, 6))
