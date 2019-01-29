#v1.0
#iterationscript 

load("variantAAtable.rda")

iteration<-function(loci){
  #subsets out OR differences smaller than 0.1 
  KDLO_List[[counter]]<-subset(KDLO_List[[counter]],  KDLO_List[[counter]][,8]>0.1)
  
  #adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
  KDLO_List[[counter]]<-unique(rbind(KDLO_List[[counter]], subset(BOLO_List[[counter]],  BOLO_List[[counter]]$Locus%in%KDLO_List[[counter]]$Locus)))[mixedorder(row.names(unique(rbind( KDLO_List[[counter]], subset( BOLO_List[[counter]],  BOLO_List[[counter]]$Locus%in% KDLO_List[[counter]]$Locus))))),]
  
  unassociated_posi<-unique(BOLO_List[[counter]]$Locus[! BOLO_List[[counter]]$Locus %in%  KDLO_List[[counter]]$Locus])
  
  if(counter==1)
  {start1<-unique( KDLO_List[[counter]]$Locus)}
  
  #for loop for all possible pair combinations 
  combinames<-sapply(start1, function(x) NULL)
  
  if(counter==1){
    for(i in 1:(length(start1)-1)){ ## range.x = 1:(N-1)
      for(j in (i+1):length(combinames)){ ## range.y = x+1:N
        if(names(combinames)[[j]]!=start1[[i]]){
          combinames[[i]][[j]]<-paste(start1[[i]],names(combinames)[[j]],sep=":")}}}
    #unlists iter0names and omits NAs to obtain all unique possible pair combinations 
    combinames<-unlist(combinames, use.names = F)[!is.na(unlist(combinames, use.names = F))]}
  
  #df for pairs -- length is number of unique pairs * 2, 
  combidf<-data.frame(variantAAtable[[loci[[1]]]][,c(1,2)], matrix("", ncol =length(rep(combinames, 2))), stringsAsFactors = F)
  
  #fills in column names 
  colnames(combidf)<-c("SampleID", "Disease", mixedsort(rep(combinames, 2)))
  
  #observes number of columns for those needed to be pasted together
  cols=c(1:length(strsplit(combinames[[1]], ":")[[1]]))
  
  #[[1]] to contain amino acid combos of TRUE/FALSE
  #[[2]] to contain amino acid combos of FALSE/TRUE
  dfAA<-sapply(1:2, function(x) NULL)
  
  #fills in element names in the lists formed in the above lists 
  for(j in 1:length(dfAA)){
    dfAA[[j]]<-sapply(combinames, function(x) NULL)}
  
  #fills in appropriate position pair combos into dfAA
  for(i in 1:length(combinames)){
    dfAA[[1]][[i]]<-apply(variantAAtable[[1]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
    dfAA[[2]][[i]]<-apply(variantAAtable[[1]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  }
  
  #fills into pair_df
  combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
  combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]

  myData<-list("combidf"=combidf, "KDLO_List"= KDLO_List[[counter]], "unassociated_posi"=unassociated_posi)
  return(myData)
  
}

myData<-NULL
myData$combidf<-variantAAtable[[1]]
stop<-FALSE
counter=0
BOLO_List <- list()

while(stop==FALSE){
  counter=counter+1
  if(counter==1){
    BOLO_List[[counter]]<-BIGDAWG(myData$combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)
    BOLO_List[[counter]]<-data.frame(lapply(as.data.frame(BOLO_List[[counter]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)
    dummy_KDLO<-as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F)[rep(seq_len(nrow(as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F))), each=nrow( BOLO_List[[counter]])),]
    dummy_KDLO[,1]<-BOLO_List[[counter]]$Locus
    dummy_KDLO[,2]<-BOLO_List[[counter]]$Allele
    for(i in 1:nrow(BOLO_List[[counter]])){
      #finds OR difference between BOLO and dummy ORs -- subs out "-", for a blank, since only evaluating absolute value of OR diff
      #adds difference to new column in BOLO 
      BOLO_List[[counter]][i,8]<-gsub("-", "", as.numeric(BOLO_List[[counter]][i,]$OR)-as.numeric(subset(subset(dummy_KDLO, grepl( BOLO_List[[counter]][i,][[1]], dummy_KDLO[,1])), grepl( BOLO_List[[counter]][i,][[2]], subset(dummy_KDLO, grepl( BOLO_List[[counter]][i,][[1]], dummy_KDLO[,1]))[,2]))[,3]))[[1]]
    }
    BOLO_List[[counter]]<-subset(BOLO_List[[counter]], BOLO_List[[counter]]$Allele!="binned")
    KDLO_List[[counter]]<-subset(BOLO_List[[counter]], BOLO_List[[counter]][,7]=="*")
  }
  
  if(counter==2){
      BOLO_List[[counter]]<-BIGDAWG(myData$combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)
      BOLO_List[[counter]]<-data.frame(lapply(as.data.frame(BOLO_List[[counter]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)
      BOLO_List[[counter]]<-subset(BOLO_List[[counter]], BOLO_List[[counter]]$Allele!="binned")
      for(i in 1:nrow(BOLO_List[[counter]])){
        BOLO_List[[counter]][i,8]<-gsub("-", "", as.numeric(BOLO_List[[counter]][i,]$OR)-as.numeric(subset(subset(myData$KDLO_List, myData$KDLO_List[,1] %in% strsplit(BOLO_List[[counter]][i,][[1]], ":")[[1]][[1]]), subset(myData$KDLO_List, myData$KDLO_List[,1] %in% strsplit(BOLO_List[[counter]][i,][[1]], ":")[[1]][[1]])$Allele %in% strsplit(BOLO_List[[counter]][i,][[2]], "~")[[1]][[1]])$OR))}
      KDLO_List[[counter]]<-subset(BOLO_List[[counter]], BOLO_List[[counter]][,7]=="*")
  }
  
  interim<-iteration("A")
  myData<-interim
  
  if(counter==2){
    print("End of motif_list analysis")
    stop=TRUE}
}

View(KDLO_List[[1]])
