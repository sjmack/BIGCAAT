##atlas formation
##version 2.0
##By:Livia Tran 

atlas<-data.frame(locus=names(framework))

atlas$UTRS <- list(c(1,17),c(1,15),c(1,17),
                   c(1,9), c(1,11), c(1,9), c(1,13), 
                   c(1,13), c(1,13), c(1,13), c(1,13))

atlas$core <- list(c(4,6),c(4,6),c(4,6),
                   4, 4, 4, 4,
                   4, 4, 4, 4)

atlas$cytoexons <- list(c(12,14,16),c(12,14),c(12,14,16),
                        8, 10, 8, c(10, 12), 
                        c(10,12), c(10,12), c(10,12), c(10,12))
atlas$leaderpep<-list(2, 2,2,
                      2,2,2,2,
                      2,2,2,2)
atlas$cdomain<-list(8,8,8,
                    8,8,8,8,
                    8,8,8,8)
atlas$introns<-list(c(3, 5, 7, 9, 11, 13, 15), c(3,5,7,9,11,13),c(3, 5, 7, 9, 11, 13, 15),
                    c(3,5,7), c(3,5,7,9), c(3,5,7), c(3,5,7,9,11),
                    c(3,5,7,9,11),c(3,5,7,9,11),c(3,5,7,9,11),c(3,5,7,9,11))
atlas$tmexons<-list(10,10,10,
                    8,8,8,8,
                    8,8,8,8)
atlas$noncore<-list(c(2, 8, 10, 12, 14, 16), c(2, 8, 10, 12), c(2, 8, 10, 12, 14, 16),
                    c(2, 6, 8), c(2, 6, 8, 10), c(2, 6, 8), c(2, 6, 8, 10, 12),
                    c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12), c(2, 6, 8, 10, 12))

setwd("~/Desktop/ltmasterscoding")
save(atlas, file="atlas.rda")
