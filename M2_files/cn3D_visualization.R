##code for obtaining full peptide sequences
#for specific motifs of interest
#to load into SWISS for protein modelling
#the protein model is downloaded as .pdb file (protein database)
#and submitted to VAST (Vector Alignment Seach Tool) to compare spatial coordinates
#of that 3D structure with the Molecular Modeling Database (MMDB)
#from there, a .CN3D file is downloaded for viewing in the CN3D program (see in 3d)
#By: Livia Tran

#load assertr library
library(assertr)

#all_gdata is a master copy of all variants of genotype data alleles and can be
#attained from the variantAA_extractor script

#peptide sequence attainment for C*03:04
#calls row 8, since that is where C*03:04 is 
iell_motif<-sapply(all_gdata$C, function(x) NULL)
for(i in 1:length(iell_motif)){
  iell_motif[[i]]<-col_concat(all_gdata[[1]][[i]][8,7:length(all_gdata[[1]][[i]])])}


#peptide sequence attainment for C*17:01
#calls row 38, since that is where C*17:01 ia 
iele_motif<-sapply(all_gdata$C, function(x) NULL)
for(i in 1:length(iele_motif)){
  iele_motif[[i]]<-col_concat(all_gdata[[1]][[i]][38,7:length(all_gdata[[1]][[i]])])}
