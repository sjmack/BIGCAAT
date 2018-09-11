#### CWDverify -- Inspecting the CWD catalogue for allele names that have been deleted. Steven J. Mack August 31, 2019 v0.1
###               This function downloads four documents from public websites (cwd.immunogenomics.org, and raw.githubusercontent.com); internet access is required.
###               The cwd200_alleles.txt, hla_nom.txt, Deleted_alleles.txt and allelelist.txt files are downloaded and processed. 
###               The function extracts the v3.0.0+ alleles that have undergone name-reassignments from the hla_nom.txt file,
###               and uses the Deleted_alleles.txt and allelelist.txt files to identify the accession numbers for the deleted and renamed alleles.
###               The accession numbers of the CWD alleles are compared to this table, and any CWD alleles that have been deleted and renamed
###               are replaced on the CWD list with the corresponding renated allele name and accession number. 

###               The function returns a two-column CWD table <accession number><allele name>.
###               Versioning information for each source file is maintained internally, and may be returned in a future version.
###               A future version should probably also verify the resolution of the remaining allele names, lengthening them as necessary. 

##                Example: newCWD <- CWDverify()

CWDverify <- function(){
require(data.table)

## Pull down the CWD catalogue
CWD <- list()
CWD$data <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",skip = 1,stringsAsFactors = FALSE,select = c(2,3))
CWD$version <- fread("https://www.uvm.edu/~igdawg/pubs/cwd200_alleles.txt",nrows = 1,stringsAsFactors = FALSE,select=1)

## Pull down the hla_nom.txt, Deleted_alleles.txt and allelelist.txt files to create a table of v3.0.0+ deleted alleles, their ACCs,their replacements, and their ACCs
deletedHLA <- list()
# Temporarily store the entire hla_nom.txt in $version
deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",skip=6, stringsAsFactors = FALSE,sep = ";", col.names = c("Locus","AlleleName","NewName","Event"),select = c(1,2,5,6))
## Exclude entries without allele name changes
deletedHLA$data <- deletedHLA$version[deletedHLA$version$NewName !="",]
# Exclude pre-db release 3.0.0 alleles
deletedHLA$data <- deletedHLA$data[grep(":",deletedHLA$data$AlleleName,fixed=TRUE),]

## Process and extract the accession numbers from the Deleted_alleles.txt file, temporarily stored in $version
deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,skip = 7,sep=",",header=TRUE,fill=TRUE)
## Below to account for one extra comma in line 106 (hopefully, can be deleted in a future release)
if(ncol(deletedHLA$version)==4) {deletedHLA$version$Description[98] <- paste(deletedHLA$version$Description[98],deletedHLA$version$V4[98],sep=" ")
deletedHLA$version <- deletedHLA$version[,1:3] }
# Store the pertinent accession numbers in the data element
deletedHLA$data$origAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$AlleleName,sep=""),deletedHLA$version$Allele)]
# Temporarily store the allelelist.txt file in $version 
deletedHLA$version <- fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",skip=6, stringsAsFactors = FALSE,sep = ",", header=TRUE)
deletedHLA$data$newAccession <- deletedHLA$version$AlleleID[match(paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep=""),deletedHLA$version$Allele)]
# overwrite the Deleted_alelles.txt files with the version information
deletedHLA$version <- cbind(fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE),fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Deleted_alleles.txt",stringsAsFactors = FALSE,nrows = 5,sep="?",header=TRUE), fread("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",nrows=5, stringsAsFactors = FALSE,sep = "?", header=TRUE))

## Match accession numbers in CWD to the Accession numbers in the deleted alleles. 
changeCWD <- match(CWD$data$`IMGT/HLA Accession Number`,deletedHLA$data$origAccession)
# Create full allele names for the new names
deletedHLA$data$NewName <- paste(deletedHLA$data$Locus,deletedHLA$data$NewName,sep="")
CWD$data[!is.na(changeCWD),] <- cbind(deletedHLA$data[changeCWD[!is.na(changeCWD)],6],deletedHLA$data[changeCWD[!is.na(changeCWD)],3])

# Rename the columns of the verified CWD table
colnames(CWD$data) <- c("Accession","AlleleName")

CWD$data
}

CWDverify()
