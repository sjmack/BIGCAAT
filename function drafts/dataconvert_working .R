dataConvert<-function(mergedcustomdata=custom_mergeddata, mapname=NULL, BIGDAWGgenotypedata, alleleListfiles, info=F){
  #requires BIGDAWG package to run loci, haplotype, and HWE tests on data
  require(BIGDAWG)
  if(info==TRUE){cat(paste("The following ‘maps’ are available in the atlas:",paste(colnames(atlas)[2:ncol(atlas)],collapse=" "),sep="\n"))}
  if(info==FALSE){
    if(any(mapname%in%colnames(atlas[,2:length(atlas)]))==FALSE){
      if(any(isTRUE(mapname=="all")==FALSE)){
        print("Error - no map name specified. Please specify a map or use 'all' to use all maps")
        mapname="all"}}
    if(any(mapname%in%colnames(atlas[,2:length(atlas)]))==TRUE){ 
      convertedlist<-list()
      #converts BIGDAWG formatted data based on feature desired, followed by immediate analysis for all
      #three available tests
      View(BDgenotypeconversion(BIGDAWGgenotypedata, alleleListfiles, custom_mergeddata[[mapname]]))
      BIGDAWG(BDgenotypeconversion(BIGDAWGgenotypedata, alleleListfiles, custom_mergeddata[[mapname]]), Run.Tests = c("HWE", "H", "L"))}
    if(any(isTRUE(mapname=="all")==TRUE)){
      convertedlist<-sapply(colnames(atlas[,2:length(atlas)]),function(x) NULL)
      #for loop for converting BIGDAWG formatted data into its GFE components based on all maps,
      #followed by immediate analysis for all tests in BIGDAWG 
      for(i in 1:length(custom_mergeddata)){
        convertedlist[[i]]<-BDgenotypeconversion(BIGDAWGgenotypedata, alleleListfiles, custom_mergeddata[[i]])}
      View(convertedlist)      
      BIGDAWG(convertedlist[[i]], HLA=F, Run.Tests = c("L", "HWE", "H"))}}}


dataConvert(mergedcustomdata = custom_mergeddata, BIGDAWGgenotypedata = HLA_data, alleleListfiles = "/Users/liviatran/Desktop/ltmasterscoding/Allelelist.3310.txt")




