#############start scripting for iteration program
require(BIGDAWG)
require(gtools)

#defines loci for BIGDAWG analysis iteration0 results 
load("variantAAtable.rda")

loci="B"

#######ITERATION 0

BOLO<-list()

BOLO[[1]]<-BIGDAWG(variantAAtable[[2]], HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

#unlists all lists to make into a dataframe 
BOLO[[1]]<-data.frame(lapply(as.data.frame(BOLO[[1]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

###only want to make the dummy for the zeroth iteration 
#makes dummy KDLO based on previous BOLO 
dummy_KDLO<-as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F)[rep(seq_len(nrow(as.data.frame(t(c("TBA-loc","TBA-allele",1.0,0.5,1.5,0.5,"NS")), stringsAsFactors = F))), each=nrow(BOLO[[1]])),]
dummy_KDLO[,1]<-BOLO[[1]]$Locus
dummy_KDLO[,2]<-BOLO[[1]]$Allele


##MAORI module 
###only want to use "dummy_KDLO" notiation for zeroth iteration 
##need to modify dummy_KDLO to just KDLO for iterations > zero 
for(i in 1:nrow(BOLO[[1]])){
  #finds OR difference between BOLO and dummy ORs -- subs out "-", for a blank, since only evaluating absolute value of OR diff
#adds difference to new column in BOLO 
  BOLO[[1]][i,8]<-gsub("-", "", as.numeric(BOLO[[1]][i,]$OR)-as.numeric(subset(subset(dummy_KDLO, grepl(BOLO[[1]][i,][[1]], dummy_KDLO[,1])), grepl(BOLO[[1]][i,][[2]], subset(dummy_KDLO, grepl(BOLO[[1]][i,][[1]], dummy_KDLO[,1]))[,2]))[,3]))[[1]]
}

KDLO<-list()

#subsets only * values for KDLO from BOLO
KDLO[[1]]<-subset(BOLO[[1]],BOLO[[1]][,7]=="*")

 #subsets out OR differences smaller than 0.1 
  KDLO[[1]]<-subset(KDLO[[1]], KDLO[[1]][,8]>0.1)

  #adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
 KDLO[[1]]<-unique(rbind(KDLO[[1]], subset(BOLO[[1]], BOLO[[1]]$Locus%in%KDLO[[1]]$Locus)))[mixedorder(row.names(unique(rbind(KDLO[[1]], subset(BOLO[[1]], BOLO[[1]]$Locus%in%KDLO[[1]]$Locus))))),]
  
 unassociated_posi<-list() 
 unassociated_posi[[1]]<-unique(BOLO[[1]]$Locus[!BOLO[[1]]$Locus %in% KDLO[[1]]$Locus])
  ####END ITERATION 0
  
  
  
  
  
##ITERATION 1   
  #start 1 should be individual amino acids from zeroth iteration
  start1<-unique(KDLO[[1]]$Locus)
  

#for loop for all possible pair combinations 
combinames<-sapply(start1, function(x) NULL)

###only want this for zeroth iteration 
for(i in 1:(length(combinames)-1)){ ## range.x = 1:(N-1)
  for(j in (i+1):length(combinames)){ ## range.y = x+1:N
    if(names(combinames)[[j]]!=start1[[i]]){
      combinames[[i]][[j]]<-paste(start1[[i]],names(combinames)[[j]],sep=":")
    }}}

#unlists iter_names and omits NAs to obtain all unique possible pair combinations 
combinames<-unlist(combinames, use.names = F)[!is.na(unlist(combinames, use.names = F))]

#df for pairs -- length is number of unique pairs * 2, 
combidf<-data.frame(variantAAtable[[2]][,c(1,2)], matrix("", ncol =length(rep(combinames,2))), stringsAsFactors = F)

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
  dfAA[[1]][[i]]<-apply(variantAAtable[[2]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  dfAA[[2]][[i]]<-apply(variantAAtable[[2]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
}

#fills into pair_df
combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]

BOLO[[2]]<-BIGDAWG(combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

BOLO[[2]]<-data.frame(lapply(as.data.frame(BOLO[[2]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

#subsets out binned alleles 
BOLO[[2]]<-subset(BOLO[[2]], (BOLO[[2]]$Allele!="binned") & (!grepl("NA", BOLO[[2]]$Allele)))

#constant[[1]][[1]] to select allele pair and first individual amino acid position 
for(i in 1:nrow(BOLO[[2]])){
  BOLO[[2]][i,8]<-gsub("-", "", as.numeric(BOLO[[2]][i,]$OR)-as.numeric(subset(subset(KDLO[[1]], KDLO[[1]][,1] %in% strsplit(BOLO[[2]][i,][[1]], ":")[[1]][[1]]), subset(KDLO[[1]], KDLO[[1]][,1] %in% strsplit(BOLO[[2]][i,][[1]], ":")[[1]][[1]])$Allele %in% strsplit(BOLO[[2]][i,][[2]], "~")[[1]][[1]])$OR))}

#subsets only * values
KDLO[[2]]<-subset(BOLO[[2]], BOLO[[2]][,7]=="*")

#subsets out OR differences smaller than 0.1 
KDLO[[2]]<-subset(KDLO[[2]], KDLO[[2]][,8]>0.1)

#adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
KDLO[[2]]<-unique(rbind(KDLO[[2]], subset(BOLO[[2]], BOLO[[2]]$Locus%in%KDLO[[2]]$Locus)))[mixedorder(row.names(unique(rbind(KDLO[[2]], subset(BOLO[[2]], BOLO[[2]]$Locus%in%KDLO[[2]]$Locus))))),]

#finds unassociated double positions 
unassociated_posi[[2]]<-unique(BOLO[[2]]$Locus[!BOLO[[2]]$Locus %in% KDLO[[2]]$Locus]) 
#some conditional for here -- is unassociated_posi the length of the original KDLO? 

##end of iteration 1





##Iteration 2
#obtains individual amino acids from splitting iter_names -- different than from just obtaining
#from KDLO$locus
start1<-unique(unlist(strsplit(combinames, ":")))

##building larger motifs for triplets 
possible_combis<-sapply(unique(KDLO[[2]]$Locus), function(x) NULL)
combinames<-NULL

#finds possible combinations by pasting names of list with singular amino acids not in that pair 
for(i in 1:length(possible_combis)){
  possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")}

#splits those triplets up and sorts them numerically to later on eliminate any duplicates 
for(j in 1:length(unlist(possible_combis))){
  combinames[[j]]<-paste(mixedsort(strsplit(unlist(possible_combis, use.names=F), ":")[[j]], decreasing=F), collapse=":")}

#eliminates duplicates
combinames<-unique(mixedsort(combinames))


#subsets out unassociated positions from sorted combinames 
for(i in 1:length(unassociated_posi[[2]])){
combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[2]][[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[2]][[i]], sep=""), combinames)))}

combidf<-data.frame(variantAAtable[[2]][,c(1,2)], matrix("", ncol =length(rep(combinames,2))), stringsAsFactors = F)

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
  dfAA[[1]][[i]]<-apply(variantAAtable[[2]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  dfAA[[2]][[i]]<-apply(variantAAtable[[2]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
}

#fills into pair_df
combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]

BOLO[[3]]<-BIGDAWG(combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

BOLO[[3]]<-data.frame(lapply(as.data.frame(BOLO[[3]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

#subsets out binned alleles 
BOLO[[3]]<-subset(BOLO[[3]], BOLO[[3]]$Allele!="binned")

###DIFFERENT BOLO STATEMENT HERE BECAUSE USING GREPL !!!!
for(i in 1:nrow(BOLO[[3]])){
  BOLO[[3]][i,8]<-gsub("-", "", as.numeric(BOLO[[3]][i,]$OR)- as.numeric(subset(subset(KDLO[[2]], KDLO[[2]][,1] %in% paste(strsplit(BOLO[[3]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[2]]$Locus, ":")[[1]]))], collapse=":")), subset(KDLO[[2]], KDLO[[2]][,1] %in% paste(strsplit(BOLO[[3]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[2]]$Locus, ":")[[1]]))], collapse=":"))$Allele %in%paste(strsplit(BOLO[[3]][i,][[2]], "~")[[1]][c(1:length(strsplit(KDLO[[2]]$Locus, ":")[[1]]))], collapse="~"))$OR))
  BOLO[[3]][i,9]<-gsub("-", "", as.numeric(BOLO[[3]][i,]$OR)-as.numeric(subset(subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[3]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[3]][i,][[1]], ":")))]]), subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[3]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[3]][i,][[1]], ":")))]])$Allele %in% strsplit(BOLO[[3]][i,][[2]], "~")[[1]][[length(unlist(strsplit(BOLO[[3]][i,][[1]], ":")))]])$OR))}

#subsets only * values
KDLO[[3]]<-subset(BOLO[[3]], BOLO[[3]][,7]=="*")

#subsets out OR differences smaller than 0.1 
KDLO[[3]]<-subset(KDLO[[3]], KDLO[[3]][,9]>0.1)
KDLO[[3]]<-subset(KDLO[[3]], KDLO[[3]][,8]>0.1)

#!!adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
KDLO[[3]]<-unique(rbind(KDLO[[3]], subset(BOLO[[3]], BOLO[[3]]$Locus%in%KDLO[[3]]$Locus)))[mixedorder(row.names(unique(rbind(KDLO[[3]], subset(BOLO[[3]], BOLO[[3]]$Locus%in%KDLO[[3]]$Locus))))),]

#finds unassociated double positions 
unassociated_posi[[3]]<-unique(BOLO[[3]]$Locus[!BOLO[[3]]$Locus %in% KDLO[[3]]$Locus]) 

##end of iteration 2 







###ITERATION 3
start1<-unique(unlist(strsplit(combinames, ":")))

##building larger motifs for triplets 
possible_combis<-sapply(unique(KDLO[[3]]$Locus), function(x) NULL)
combinames<-NULL

#finds possible combinations by pasting names of list with singular amino acids not in that pair 
for(i in 1:length(possible_combis)){
  possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")}

#splits those triplets up and sorts them numerically to later on eliminate any duplicates 
for(j in 1:length(unlist(possible_combis))){
  combinames[[j]]<-paste(mixedsort(strsplit(unlist(possible_combis, use.names=F), ":")[[j]], decreasing=F), collapse=":")}

#eliminates duplicates
combinames<-unique(mixedsort(combinames))

View(combinames)

#subsets out unassociated positions from sorted combinames 
for(i in 1:length(unassociated_posi[[3]])){
  combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[3]][[i]], sep=""), combinames)))}

for(j in 1:length(unassociated_posi[[2]])){
  combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[2]][[j]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[2]][[j]], sep=""), combinames)))}

combidf<-data.frame(variantAAtable[[2]][,c(1,2)], matrix("", ncol =length(rep(combinames, 2))), stringsAsFactors = F)

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
  dfAA[[1]][[i]]<-apply(variantAAtable[[2]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  dfAA[[2]][[i]]<-apply(variantAAtable[[2]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
}

#fills into pair_df
combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]

BOLO[[4]]<-BIGDAWG(combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

BOLO[[4]]<-data.frame(lapply(as.data.frame(BOLO[[4]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

#subsets out binned alleles 
BOLO[[4]]<-subset(BOLO[[4]], BOLO[[4]]$Allele!="binned")


#constant[[1]][[1]] to select allele pair and first individual amino acid position 
for(i in 1:nrow(BOLO[[4]])){
    BOLO[[4]][i,8]<-gsub("-", "", as.numeric(BOLO[[4]][i,]$OR)- as.numeric(subset(subset(KDLO[[3]], KDLO[[3]][,1] %in% paste(strsplit(BOLO[[4]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[3]]$Locus, ":")[[1]]))], collapse=":")), subset(KDLO[[3]], KDLO[[3]][,1] %in% paste(strsplit(BOLO[[4]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[3]]$Locus, ":")[[1]]))], collapse=":"))$Allele %in%paste(strsplit(BOLO[[4]][i,][[2]], "~")[[1]][c(1:length(strsplit(KDLO[[3]]$Locus, ":")[[1]]))], collapse="~"))$OR))
    BOLO[[4]][i,9]<-gsub("-", "", as.numeric(BOLO[[4]][i,]$OR)-as.numeric(subset(subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[4]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[4]][i,][[1]], ":")))]]), subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[4]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[4]][i,][[1]], ":")))]])$Allele %in% strsplit(BOLO[[4]][i,][[2]], "~")[[1]][[length(unlist(strsplit(BOLO[[4]][i,][[1]], ":")))]])$OR))}


#subsets only * values
KDLO[[4]]<-subset(BOLO[[4]], BOLO[[4]][,7]=="*")

#subsets out OR differences smaller than 0.1 
KDLO[[4]]<-subset(KDLO[[4]], KDLO[[4]][,8]>0.1)

KDLO[[4]]<-subset(KDLO[[4]], KDLO[[4]][,9]>0.1)

#adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
KDLO[[4]]<-unique(rbind(KDLO[[4]], subset(BOLO[[4]], BOLO[[4]]$Locus%in%KDLO[[4]]$Locus)))[mixedorder(row.names(unique(rbind(KDLO[[4]], subset(BOLO[[4]], BOLO[[4]]$Locus%in%KDLO[[4]]$Locus))))),]

#finds unassociated double positions 
unassociated_posi[[4]]<-unique(BOLO[[4]]$Locus[!BOLO[[4]]$Locus %in% KDLO[[4]]$Locus]) 


View(unassociated_posi[[3]])
View(combinames)

#iteration 5

##building larger motifs for triplets 
possible_combis<-sapply(unique(KDLO[[4]]$Locus), function(x) NULL)
combinames<-NULL

#finds possible combinations by pasting names of list with singular amino acids not in that pair 
for(i in 1:length(possible_combis)){
  possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")}

#splits those triplets up and sorts them numerically to later on eliminate any duplicates 
for(j in 1:length(unlist(possible_combis))){
  combinames[[j]]<-paste(mixedsort(strsplit(unlist(possible_combis, use.names=F), ":")[[j]], decreasing=F), collapse=":")}

#eliminates duplicates
combinames<-unique(mixedsort(combinames))

#subsets out unassociated positions from sorted combinames 
for(i in 1:length(unassociated_posi[[4]])){
  combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[4]][[i]], sep=""), combinames))}

for(i in 1:length(unassociated_posi[[3]])){
  combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[3]][[i]], sep=""), combinames))}

for(i in 1:length(unassociated_posi[[2]])){
  combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[2]][[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[2]][[i]], sep=""), combinames)))}


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
  dfAA[[1]][[i]]<-apply(variantAAtable[[2]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  dfAA[[2]][[i]]<-apply(variantAAtable[[2]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
}

#fills into pair_df
combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]

BOLO[[5]]<-BIGDAWG(combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

BOLO[[5]]<-data.frame(lapply(as.data.frame(BOLO[[5]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

#subsets out binned alleles 
BOLO[[5]]<-subset(BOLO[[5]], BOLO[[5]]$Allele!="binned")

###DIFFERENT BOLO STATEMENT HERE BECAUSE USING GREPL !!!!
for(i in 1:nrow(BOLO[[5]])){
  BOLO[[5]][i,8]<-gsub("-", "", as.numeric(BOLO[[5]][i,]$OR)- as.numeric(subset(subset(KDLO[[4]], KDLO[[4]][,1] %in% paste(strsplit(BOLO[[5]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[4]]$Locus, ":")[[1]]))], collapse=":")), subset(KDLO[[4]], KDLO[[4]][,1] %in% paste(strsplit(BOLO[[5]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[4]]$Locus, ":")[[1]]))], collapse=":"))$Allele %in%paste(strsplit(BOLO[[5]][i,][[2]], "~")[[1]][c(1:length(strsplit(KDLO[[4]]$Locus, ":")[[1]]))], collapse="~"))$OR))
  BOLO[[5]][i,9]<-gsub("-", "", as.numeric(BOLO[[5]][i,]$OR)-as.numeric(subset(subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[5]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[5]][i,][[1]], ":")))]]), subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[5]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[5]][i,][[1]], ":")))]])$Allele %in% strsplit(BOLO[[5]][i,][[2]], "~")[[1]][[length(unlist(strsplit(BOLO[[5]][i,][[1]], ":")))]])$OR))}

#subsets only * values
KDLO[[5]]<-subset(BOLO[[5]], BOLO[[5]][,7]=="*")

#subsets out OR differences smaller than 0.1 
KDLO[[5]]<-subset(KDLO[[5]], KDLO[[5]][,9]>0.1)
KDLO[[5]]<-subset(KDLO[[5]], KDLO[[5]][,8]>0.1)

#!!adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
KDLO[[5]]<-unique(rbind(KDLO[[5]], subset(BOLO[[5]], BOLO[[5]]$Locus%in%KDLO[[5]]$Locus)))[mixedorder(row.names(unique(rbind(KDLO[[5]], subset(BOLO[[5]], BOLO[[5]]$Locus%in%KDLO[[5]]$Locus))))),]

#finds unassociated double positions 
unassociated_posi[[5]]<-unique(BOLO[[5]]$Locus[!BOLO[[5]]$Locus %in% KDLO[[5]]$Locus]) 



##ITERATION 6
start1<-unique(unlist(strsplit(combinames, ":")))

##building larger motifs for triplets 
possible_combis<-sapply(unique(KDLO[[5]]$Locus), function(x) NULL)
combinames<-NULL

#finds possible combinations by pasting names of list with singular amino acids not in that pair 
for(i in 1:length(possible_combis)){
  possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")}

#splits those triplets up and sorts them numerically to later on eliminate any duplicates 
for(j in 1:length(unlist(possible_combis))){
  combinames[[j]]<-paste(mixedsort(strsplit(unlist(possible_combis, use.names=F), ":")[[j]], decreasing=F), collapse=":")}

#eliminates duplicates
combinames<-unique(mixedsort(combinames))

#subsets out unassociated positions from sorted combinames 
for(i in 1:length(unassociated_posi[[5]])){
  combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[5]][[i]], sep=""), combinames))}

for(i in 1:length(unassociated_posi[[4]])){
  combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[4]][[i]], sep=""), combinames))}

for(i in 1:length(unassociated_posi[[3]])){
  combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[3]][[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[3]][[i]], sep=""), combinames)))}

for(i in 1:length(unassociated_posi[[2]])){
  combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[2]][[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[2]][[i]], sep=""), combinames)))}

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
  dfAA[[1]][[i]]<-apply(variantAAtable[[2]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  dfAA[[2]][[i]]<-apply(variantAAtable[[2]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
}

#fills into pair_df
combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]

BOLO[[6]]<-BIGDAWG(combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

BOLO[[6]]<-data.frame(lapply(as.data.frame(BOLO[[6]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

#subsets out binned alleles 
BOLO[[6]]<-subset(BOLO[[6]], BOLO[[6]]$Allele!="binned")


#constant[[1]][[1]] to select allele pair and first individual amino acid position 

for(i in 1:nrow(BOLO[[6]])){
  BOLO[[6]][i,8]<-gsub("-", "", as.numeric(BOLO[[6]][i,]$OR)- as.numeric(subset(subset(KDLO[[5]], KDLO[[5]][,1] %in% paste(strsplit(BOLO[[6]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[5]]$Locus, ":")[[1]]))], collapse=":")), subset(KDLO[[5]], KDLO[[5]][,1] %in% paste(strsplit(BOLO[[6]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[5]]$Locus, ":")[[1]]))], collapse=":"))$Allele %in%paste(strsplit(BOLO[[6]][i,][[2]], "~")[[1]][c(1:length(strsplit(KDLO[[5]]$Locus, ":")[[1]]))], collapse="~"))$OR))
  BOLO[[6]][i,9]<-gsub("-", "", as.numeric(BOLO[[6]][i,]$OR)-as.numeric(subset(subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[6]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[6]][i,][[1]], ":")))]]), subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[6]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[6]][i,][[1]], ":")))]])$Allele %in% strsplit(BOLO[[6]][i,][[2]], "~")[[1]][[length(unlist(strsplit(BOLO[[6]][i,][[1]], ":")))]])$OR))}


#subsets only * values
KDLO[[6]]<-subset(BOLO[[6]], BOLO[[6]][,7]=="*")

#subsets out OR differences smaller than 0.1 
KDLO[[6]]<-subset(KDLO[[6]], KDLO[[6]][,8]>0.1)

KDLO[[6]]<-subset(KDLO[[6]], KDLO[[6]][,9]>0.1)


#adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
KDLO[[6]]<-unique(rbind(KDLO[[6]], subset(BOLO[[6]], BOLO[[6]]$Locus%in%KDLO[[6]]$Locus)))[mixedorder(row.names(unique(rbind(KDLO[[6]], subset(BOLO[[6]], BOLO[[6]]$Locus%in%KDLO[[6]]$Locus))))),]

#finds unassociated double positions 
unassociated_posi[[6]]<-unique(BOLO[[6]]$Locus[!BOLO[[6]]$Locus %in% KDLO[[6]]$Locus]) 


##iteration 7 
start1<-unique(unlist(strsplit(combinames, ":")))

##building larger motifs for triplets 
possible_combis<-sapply(unique(KDLO[[6]]$Locus), function(x) NULL)
combinames<-NULL

#finds possible combinations by pasting names of list with singular amino acids not in that pair 
for(i in 1:length(possible_combis)){
  possible_combis[[i]]<-paste(names(possible_combis[i]), unique(start1[which(start1%in%strsplit(names(possible_combis[i]), ":")[[1]]==FALSE)]), sep=":")}

#splits those triplets up and sorts them numerically to later on eliminate any duplicates 
for(j in 1:length(unlist(possible_combis))){
  combinames[[j]]<-paste(mixedsort(strsplit(unlist(possible_combis, use.names=F), ":")[[j]], decreasing=F), collapse=":")}

#eliminates duplicates
combinames<-unique(mixedsort(combinames))

#subsets out unassociated positions from sorted combinames 
if(length(unassociated_posi[[6]])!=0){
for(i in 1:length(unassociated_posi[[6]])){
  combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[6]][[i]], sep=""), combinames))}
}

for(i in 1:length(unassociated_posi[[5]])){
  combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[5]][[i]], sep=""), combinames))}

for(i in 1:length(unassociated_posi[[4]])){
  combinames<-subset(combinames, !grepl(paste("^", unassociated_posi[[4]][[i]], sep=""), combinames))}

for(i in 1:length(unassociated_posi[[3]])){
  combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[3]][[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[3]][[i]], sep=""), combinames)))}

for(i in 1:length(unassociated_posi[[2]])){
  combinames<-subset(combinames, (!grepl(paste("^", unassociated_posi[[2]][[i]], sep=""), combinames)) & (!grepl(paste(":", unassociated_posi[[2]][[i]], sep=""), combinames)))}

combidf<-data.frame(variantAAtable[[2]][,c(1,2)], matrix("", ncol =length(rep(combinames, 2))), stringsAsFactors = F)

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
  dfAA[[1]][[i]]<-apply(variantAAtable[[2]][c(TRUE, FALSE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
  dfAA[[2]][[i]]<-apply(variantAAtable[[2]][c(FALSE, TRUE)][strsplit(combinames, ":")[[i]]][,cols], 1, paste, collapse = "~")
}

#fills into pair_df
combidf[,3:length(combidf)][,c(TRUE,FALSE)]<-dfAA[[1]]
combidf[,3:length(combidf)][,c(FALSE,TRUE)]<-dfAA[[2]]

BOLO[[7]]<-BIGDAWG(combidf, HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)

BOLO[[7]]<-data.frame(lapply(as.data.frame(BOLO[[7]]$L$Set1$OR), function(x) unlist(x)), stringsAsFactors = F)

#subsets out binned alleles 
BOLO[[7]]<-subset(BOLO[[7]], BOLO[[7]]$Allele!="binned")

#constant[[1]][[1]] to select allele pair and first individual amino acid position 

for(i in 1:nrow(BOLO[[7]])){
  BOLO[[7]][i,8]<-gsub("-", "", as.numeric(BOLO[[7]][i,]$OR)- as.numeric(subset(subset(KDLO[[6]], KDLO[[6]][,1] %in% paste(strsplit(BOLO[[7]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[6]]$Locus, ":")[[1]]))], collapse=":")), subset(KDLO[[6]], KDLO[[6]][,1] %in% paste(strsplit(BOLO[[7]][i,][[1]], ":")[[1]][c(1:length(strsplit(KDLO[[6]]$Locus, ":")[[1]]))], collapse=":"))$Allele %in%paste(strsplit(BOLO[[7]][i,][[2]], "~")[[1]][c(1:length(strsplit(KDLO[[6]]$Locus, ":")[[1]]))], collapse="~"))$OR))
  BOLO[[7]][i,9]<-gsub("-", "", as.numeric(BOLO[[7]][i,]$OR)-as.numeric(subset(subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[7]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[7]][i,][[1]], ":")))]]), subset(KDLO[[1]], KDLO[[1]]$Locus %in% strsplit(BOLO[[7]][i,][[1]], ":")[[1]][[length(unlist(strsplit(BOLO[[7]][i,][[1]], ":")))]])$Allele %in% strsplit(BOLO[[7]][i,][[2]], "~")[[1]][[length(unlist(strsplit(BOLO[[7]][i,][[1]], ":")))]])$OR))}

#subsets only * values
KDLO[[6]]<-subset(BOLO[[6]], BOLO[[6]][,7]=="*")

#subsets out OR differences smaller than 0.1 
KDLO[[6]]<-subset(KDLO[[6]], KDLO[[6]][,8]>0.1)

KDLO[[6]]<-subset(KDLO[[6]], KDLO[[6]][,9]>0.1)
View(KDLO[[6]])

#adds in positions from original BOLO that were previously eliminated because of NS or <0.1 variant  
KDLO[[6]]<-unique(rbind(KDLO[[6]], subset(BOLO[[6]], BOLO[[6]]$Locus%in%KDLO[[6]]$Locus)))[mixedorder(row.names(unique(rbind(KDLO[[6]], subset(BOLO[[6]], BOLO[[6]]$Locus%in%KDLO[[6]]$Locus))))),]

#finds unassociated double positions 
unassociated_posi[[6]]<-unique(BOLO[[6]]$Locus[!BOLO[[6]]$Locus %in% KDLO[[6]]$Locus]) 


interim$combinames

paste(":", UMLO_list[[counter-2]][[j]], sep="")


subset(interim$combinames, !grepl(paste(":", UMLO_list[[counter-2]][[1]], sep=""), interim$combinames))




