GenotypeGFEmatch <- function(allelefiles,cwddata,gfepath, genotypedatafiles){
hlamerged<-CWDrestriction(allelefiles,cwddata,gfepath)
hlafields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(hlamerged[,1],"\\*")), nrow=nrow(hlamerged), byrow=T), stringsAsFactors=FALSE)) [2,],":")
genodata<-filemerge(genotypedatafiles, "genodata" , 0, 0)
genofields<-strsplit(t(as.data.frame(matrix(unlist(strsplit(genodata[,1],"\\*")), nrow=nrow(genodata), byrow=T), stringsAsFactors=FALSE)) [2,],":")
genodata$gfeCWD<-ifelse(hlamerged$CWD[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="CWD", 
                          hlamerged$gfe[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                          )], ifelse(hlamerged$CWD[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*"))]=="NON-CWD", 
                                     hlamerged$gfe[match(paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", genodata$genodata), paste(sapply(genofields, "[", 1), sapply(genofields, "[", 2), sep=":"), sep="*"), paste(gsub(".*[HLA-]([^*]+)[*].*", "\\1", hlamerged$allelename), paste(sapply(hlafields, "[", 1), sapply(hlafields, "[", 2), sep=":"), sep="*")
                                     )], "")) 
  return(genodata)}

GenotypeGFEmatch("Allelelist.3310.txt","cwd200_alleles.txt","/Users/liviatran/Desktop/ltmasterscoding/HLA","/Users/liviatran/Desktop/ltmasterscoding/test")
