#v 1.0

loci="A"

#######ITERATION 0 
motif_list<-c(0,2)

load("variantAAtable.rda")

iteration<-function(loci, myData, motif_list, counter){
BOLO<-sapply(loci, function(x) NULL)

BOLO<-BIGDAWG(myData, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

#unlists all lists in the dataframe
BOLO<-data.frame(lapply(as.data.frame(BOLO$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

if(counter==0){
  #makes dummy KDLO based on previous BOLO 
  dummy_KDLO<-as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F)[rep(seq_len(nrow(as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F))), each=nrow(BOLO)),]
  dummy_KDLO[,1]<-BOLO$Locus
  dummy_KDLO[,2]<-BOLO$Allele

##MAORI module 
for(i in 1:nrow(BOLO)){
  #finds OR difference between BOLO and dummy ORs -- subs out "-", for a blank, since only evaluating absolute value of OR diff
  #adds difference to new column in BOLO 
  BOLO[i,8]<-gsub("-", "", as.numeric(BOLO[i,]$OR)-as.numeric(subset(subset(dummy_KDLO, grepl(BOLO[i,][[1]], dummy_KDLO[,1])), grepl(BOLO[i,][[2]], subset(dummy_KDLO, grepl(BOLO[i,][[1]], dummy_KDLO[,1]))[,2]))[,3]))[[1]]
}}

if(counter==1){
 for(i in 1:nrow(BOLO)){
  BOLO[i,8]<-gsub("-", "", as.numeric(BOLO[i,]$OR)-as.numeric(subset(subset(KDLO, KDLO[,1] %in% strsplit(BOLO[i,][[1]], ":")[[1]][[1]]), subset(KDLO, KDLO[,1] %in% strsplit(BOLO[i,][[1]], ":")[[1]][[1]])$Allele %in% strsplit(BOLO[i,][[2]], "~")[[1]][[1]])$OR))}
}}


#subsets only * values for KDLO from BOLO
KDLO<-subset(BOLO,BOLO[,7]=="*")

#subsets out OR differences smaller than 0.1 
KDLO<-subset(KDLO, KDLO[,8]>0.1)

#adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
KDLO<-unique(rbind(KDLO, subset(BOLO, BOLO$Locus%in%KDLO$Locus)))[mixedorder(row.names(unique(rbind(KDLO, subset(BOLO, BOLO$Locus%in%KDLO$Locus))))),]

unassociated_posi<-unique(BOLO$Locus[!BOLO$Locus %in% KDLO$Locus])

if(counter==0)
  {start1<-unique(KDLO$Locus)}

if(counter==1){
start1<-"1"}

#for loop for all possible pair combinations 
combinames<-sapply(start1, function(x) NULL)

if(counter==0){
  for(i in 1:(length(start1)-1)){ ## range.x = 1:(N-1)
  for(j in (i+1):length(combinames)){ ## range.y = x+1:N
  if(names(combinames)[[j]]!=start1[[i]]){
  combinames[[i]][[j]]<-paste(start1[[i]],names(combinames)[[j]],sep=":")}}}
  #unlists iter0names and omits NAs to obtain all unique possible pair combinations 
  combinames<-unlist(combinames, use.names = F)[!is.na(unlist(combinames, use.names = F))]}

  if(counter==1){
  possible_combis<-sapply(unique(KDLO$Locus), function(x) NULL)
  
  #finds possible combinations by pasting names of list with singular amino acids not in that pair 
  #for(i in 1:length(possible_combis)){
  possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")
  
  #splits those triplets up and sorts them numerically to later on eliminate any duplicates 
  #for(j in 1:length(unlist(possible_combis))){
   combinames[[j]]<-paste(mixedsort(strsplit(unlist(possible_combis, use.names=F), ":")[[j]], decreasing=F), collapse=":")}


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



myData<-list("KDLO"=KDLO, "BOLO"=BOLO, "combidf"=combidf, "unassociated_posi"=unassociated_posi, "combinames"=combinames)
return(myData)


}


myData$combidf<-variantAAtable[[1]]
stop<-FALSE

counter=0
while(stop==FALSE){
  interim<-iteration("A", myData$combidf, c(0,2), 0)
  counter=counter+1

  myData$combidf<-interim$combidf
  myData<-interim
  

  if(counter==2){
    print("End of motif_list analysis")
    stop=TRUE}
}





if(any(motif_list==2)){
  start1<-unique(unlist(strsplit(combinames, ":")))
  possible_combis<-sapply(unique(combinames), function(x) NULL)
  for(i in 1:length(possible_combis)){
    possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")}
  possible_combis<-strsplit(unlist(possible_combis, use.names=F), ":")
  for(k in 1:length(possible_combis)){
    possible_combis[[k]]<-paste(mixedsort(possible_combis[[k]]),collapse=":")}
  combinames<-mixedsort(unique(unlist(possible_combis, use.names=F)))
  #subsets out unassociated positions from sorted combinames 
  for(i in 1:length(unassociated_posi)){
    combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[i]], sep=""), combinames))}
  unassociated_posi<-unique(BOLO$Locus[!BOLO$Locus %in% KDLO$Locus])
  for(i in 1:length(unassociated_posi)){
    combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[i]], sep=""), combinames))}
}

if(motif_list[[2]]==2){
  for(i in 1:nrow(BOLO)){
    #finds OR difference between BOLO and dummy ORs -- subs out "-", for a blank, since only evaluating absolute value of OR diff
    #adds difference to new column in BOLO 
    BOLO[i,8]<-gsub("-", "", as.numeric(BOLO[i,]$OR)-as.numeric(subset(subset(KDLO, grepl(BOLO[i,][[1]], KDLO[,1])), grepl(BOLO[i,][[2]], subset(KDLO, grepl(BOLO[i,][[1]], KDLO[,1]))[,2]))[,3]))[[1]]
  }}



