#Correspondence table 
#for determining differences between actual and alignment starts and ends of amino acid sequences
#v 0.1 
#By: Livia Tran 


corr_table<-list()
loci<-c("A", "B", "C")


corr_table[[loci[[1]]]]<- data.frame("correspondence"=T,
                                     "actual_start"=-24,
                                     "actual_end"=365,
                                     "alignment_start"=1,
                                     "alignment_end"=341, stringsAsFactors = FALSE)


corr_table[[loci[[2]]]]<- data.frame("correspondence"=c(T, F, T),
                                     "actual_start"=c(-23,1,295),
                                     "actual_end"=c(294,1,349), 
                                     "alignment_start"=c(1, 295, 296),
                                     "alignment_end"=c(294, 295, 350), stringsAsFactors = FALSE)



corr_table[[loci[[3]]]]<- data.frame("correspondence"=c(T,F,T), 
                                     "actual_start"=c(-24,1,301),
                                     "actual_end"=c(300,6, 355),
                                     "alignment_start"=c(1, 325, 331),
                                     "alignment_end"=c(324, 330, 405), stringsAsFactors = FALSE)

View(corr_table[[3]])
