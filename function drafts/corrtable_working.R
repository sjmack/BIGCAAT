#correspondence_table working script
#October 30th, 2018
#v 0.1


#Note: Currently in working progress for HLA-C only
loci="C"

#requires stringr package 
require(stringr)

##script to extract the number of characters in the reference sequence for later usage
#in determining the number of characters in the alignment sequence, including any inDels
start<-end<-list()
downloaded<-list()
for(i in 1:length(loci)) {
  downloaded[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(loci,"_prot.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
  downloaded[[loci[i]]] <- head(downloaded[[i]],-3)
  downloaded[[loci[i]]] <- tail(downloaded[[i]],-7)
  downloaded[[loci[i]]]<-sub(" ", "", downloaded[[i]])
  downloaded[[loci[i]]] <-paste(substr(downloaded[[i]],1,regexpr(" ",text = downloaded[[i]],fixed = TRUE)), gsub(" ","",substr(downloaded[[i]],regexpr(" ",text = downloaded[[i]],fixed = TRUE),nchar(downloaded[[i]]))),sep = "")
  downloaded[[loci[i]]]  <- strsplit(downloaded[[i]]," ", fixed=T)
  downloaded[[loci[i]]] <- do.call(rbind,downloaded[[i]])
  colnames(downloaded[[i]])<-c(paste(loci[[i]], "alleles", sep=""), "pepseq")
  start[[loci[i]]]<-as.numeric(which(downloaded[[i]][,1]=="Prot"))
  end[[loci[i]]] <- as.numeric(c(start[[i]][2:length(start[[i]])]-1,nrow(downloaded[[i]])))
  for(j in 1:length(start[[i]])){
    if(nrow(downloaded[[i]][start[[i]][j]:end[[i]][j],])!=nrow(downloaded[[i]][start[[i]][1]:end[[i]][1],]))
    {x<-as.data.frame(downloaded[[i]][,1][start[[i]][1]:end[[i]][1]][-c(1,2)], stringsAsFactors = F)
    colnames(x)<-paste(loci[[i]], "alleles", sep="")
    x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(downloaded[[i]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
    y<-data.frame(tail(downloaded[[i]],1), stringsAsFactors = F)
    x$pepseq[match(y[,1], x[,1])]<-y$pepseq
    downloaded[[i]]<-as.matrix(rbind(head(downloaded[[i]], -1), x))
    start[[loci[i]]]<-as.numeric(which(downloaded[[i]][,1]=="Prot"))
    end[[loci[i]]] <- as.numeric(c(start[[i]][2:length(start[[i]])]-1,nrow(downloaded[[i]])))}}}

downloaded_segments<-sapply(loci, function(x) NULL)
for(i in 1:length(loci)){
  for(j in 1:length(start[[i]])){
    downloaded_segments[[i]]<-cbind(downloaded_segments[[i]], downloaded[[i]][start[[i]][j]:end[[i]][j],])}
  cols<-seq(0, ncol(downloaded_segments[[i]]), by=2)
  downloaded_segments[[i]]<-cbind(downloaded_segments[[i]][,1], apply(downloaded_segments[[i]][,cols], 1 ,paste, collapse = ""))
  rownames(downloaded_segments[[i]])=NULL}


#downloads text files from ANHIG -- preserves white spaces and structure
#removes non-pertinent header and footer detail 
txtfile <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/",paste(loci,"_prot.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
txtfile <- tail(txtfile,-7)
txtfile<-head(txtfile, -3)

#difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing 
diff<-(countSpaces(txtfile[2])[1])-(countSpaces(txtfile[1])[2])+1

#determines positions of "Prot"
prot_start<-as.numeric(grep("Prot", txtfile))
#prot_end <- c(prot_coord[2:length(prot_coord)]-1,nrow(txtfile))

#counts number of characters in the very last allele to add onto the last Prot enumeration block
#obtains end length 
end_char<-nchar(sapply(strsplit(gsub(" ", "", sub(" ", "~", str_squish(tail(txtfile, 1)))), "~"), "[", 2))-1

prot_extractions<- refblock_number<-NULL

#extracts rows with Prot and reference sequence position information 
#extracts only relevant reference sequence positions
#NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
#as the actual sequence start will always be 1 
for (i in 1:length(prot_start)){
prot_extractions[[i]]<-txtfile[prot_start[[i]]]
prot_extractions[[i]]<-str_squish(prot_extractions[[i]])
prot_extractions[[i]]<-strsplit(prot_extractions[[i]], " ")
refblock_number[[i]]<-sapply(prot_extractions[[i]], "[", 2)
}

#determines the alignment start by adding -30 to the difference between white spaces found above 
alignment_start<-as.numeric(refblock_number[[1]])+diff

#creates a matrix where number of columns is based on the number of characters in the alignment reference sequence
#including inDels 
df<-matrix(, nrow = 2, ncol = as.numeric(nchar(downloaded_segments[[1]][,2][3])))

#inputs numeric 1:the number of characters in the reference sequence for the **actual sequence 
df[1,]<-(1:as.numeric(nchar(downloaded_segments[[1]][,2][3])))

#determines alignment length based on the last Prot enumeration + the end character 
alignment_length<-as.numeric(tail(refblock_number, 1))+end_char

#pastes -24 to 375 together in sequential order, with inDels accounted for 
#captures output as "w"
w<-capture.output(cat(alignment_start:(as.numeric(refblock_number[4])+24), paste("inDel", seq(1:6), sep=""), (as.numeric(refblock_number[4])+24+1):alignment_length))

#splits string formed by cat for separate character variables
alignment_seq<-as.character(unlist(strsplit(w, " ")))

#eliminates "0", as the alignment sequence from ANHIG does not have 0
alignment_seq<-alignment_seq[-which(alignment_seq == 0)]

#inserts into matrix previously formed 
df[2,]<-alignment_seq




##function to count spaces in between regions of interest
#to determine where  start for the alignment sequence 
countSpaces <- function(x){
  counter <- 0
  coll <- numeric()
  vec <- strsplit(x," ")[[1]]
  for(i in 1:length(vec)){
    if (vec[i]==""){
      counter <- counter+1
    }
    else{
      if (counter!=0) coll <- c(coll,counter)
      counter <- 1
    }
  }
  coll
}
