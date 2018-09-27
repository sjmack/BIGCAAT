dataConvert<-function(mergedcustomdata=custom_mergeddata, mapname="all", BIGDAWGgenotypedata, alleleListfiles){
  if(any(mapname%in%colnames(atlas[,2:length(atlas)]))==TRUE){
  start=1
  end=length(custom_mergeddata[[mapname]])
  convertedlist<-list()
  i=mapname}
      if(any(mapname=="all")){
        convertedlist<-sapply(colnames(atlas[,2:length(atlas)]),function(x) NULL)
        start=1
        end=ncol(atlas[,2:length(atlas)])}
       for (i in start:end){
        convertedlist[[i]]<-BDgenotypeconversion(HLA_data, "/Users/liviatran/Desktop/ltmasterscoding/Allelelist.3310.txt", custom_mergeddata[[i]])}
  return(convertedlist)}


dataConvert(custom_mergeddata, "fiveUTR", HLA_data, "/Users/liviatran/Desktop/ltmasterscoding/Allelelist.3310.txt")




