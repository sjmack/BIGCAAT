#defines loci for BIGDAWG analysis iteration0 results 

loci="A"

iterations<-sapply(loci, function(x) NULL)

for(i in 1:length(loci)){
  #runs BIGDAWG on all loci 
  iterations[[loci[[i]]]]<-BIGDAWG(variantAAtable[[i]], HLA=F, Run.Tests="L", Missing = "ignore", Return=T, Output = F)
  
  #subsets iteration 0 variantAAtable amino acids to those that were significant in the first round of BIGDAWG
  #analysis 
  for(j in 1:length(iterations[[loci[[i]]]])){
    iterations[[loci[i]]]<-as.data.frame(subset(iterations[[j]]$L$Set1$OR, iterations[[j]]$L$Set1$OR[,7]=="*"), stringsAsFactors=FALSE)}
  
  iteration_names<-sapply(iterations[[loci[[i]]]]$Locus, function(x) NULL)
  
  setStart<-iterations[[loci[[i]]]]$Locus
  
  for(d in 1:(length(setStart)-1)) { ## range.x = 1:(N-1)
    for(b in (d+1):length(setStart)){ ## range.y = x+1:N
      if(setStart[[b]]!=setStart[[d]]){
        iteration_names[[b]][[d]]<-paste(setStart[[d]],setStart[[b]],sep=":")
      }
    }}
}
  
  iterations_df<-data.frame(variantAAtable[[loci[[1]]]][,c(1,2)], matrix("", ncol =length(mixedsort(rep(unique(unlist(pairNames)),2)))), stringsAsFactors = F)
  colnames(iterations_df)<-c("SampleID", "Disease", mixedsort(rep(unique(unlist(iteration_names)),2)))

##these amino acid pairs move on   
  pairPass<-sapply(1:(nrow(iterations[[1]])-1), function(x) NULL)
  
  names(pairPass)<-iterations[[1]]$Locus[1:(nrow(iterations[[1]])-1)]
  
  for(j in 1:(length(pairPass))){
    pairPass[[j]]<-subset(subset(iteration1BDresults$L$Set1$OR, grepl(paste("^", paste(iterations[[1]][j,][[1]], ":", sep=""), sep=""), iteration1BDresults$L$Set1$OR[,1])), grepl(paste("^", paste(iterations[[1]][j,][[2]], sep="~"), sep=""), subset(iteration1BDresults$L$Set1$OR, grepl(paste("^", paste(iterations[[1]][j,][[1]], ":", sep=""), sep=""), iteration1BDresults$L$Set1$OR[,1]))[,2]))
    pairPass[[j]]<-cbind(pairPass[[j]], "ORdiff"=gsub("-", "", as.numeric(iterations[[1]][j,][[3]])-as.numeric(pairPass[[j]][,3]))) 
    pairPass[[j]]<-subset(pairPass[[j]], pairPass[[j]][,8]>0.1)
  }
  pairPass<-pairPass[sapply(pairPass, nrow)>0]  
  
  
  
  
  
  ###
  
  thing=NA
  for(i in 1:(length(setStart)-1)) { ## range.x = 1:(N-1)
    for(j in (i+1):length(setStart)){ ## range.y = x+1:N
      if(any(is.na(thing)==T)){
        if(setStart[[j]]!=setStart[[i]]){
          print(paste(setStart[[i]],setStart[[j]],sep=":"))}
      }
    }
  }
  
  #modified algorithm for pairs
  thing=0
  next_iter=NULL
  meep<-sapply(mixedsort(unique(unlist(pairNames))), function(x) NULL)
  if(any(is.na(thing)==F)){
    for(i in 1:length(meep)){
      meep[[i]]<-paste(names(meep[i]), unique(setStart[which(setStart%in%strsplit(names(meep[i]), ":")[[1]]==FALSE)]), sep=":")}
    for(i in 1:length(strsplit(unlist(meep, use.names=F), ":"))){
      next_iter[[i]]<-paste(mixedsort(strsplit(unlist(meep, use.names=F), ":")[[i]], decreasing=F), collapse=":")}
    next_iter<-unique(mixedsort(next_iter))
  }