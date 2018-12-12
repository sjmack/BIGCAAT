#variant match code 



#creates a new list with 2 columns, since each locus in gdata has 2 columns
#variant_match<-sapply(1:2, function(x) NULL)

#further adds on to the previously created variable, where the final product is 
#a list with 2 elements
#where each element contains a list of 831 elements, i.e. the number of rows in gdata
#where each of the 831 elements contains a list of 145 elements, i.e. the length of position_parsed
#for(i in 1:length(variant_match)){
# variant_match[[i]]<-sapply(nrow(gdata), function(x) NULL)
#for(j in 1:nrow(gdata)){
# variant_match[[i]][[j]]<-sapply(position_parsed, function(x) NULL)
#}
#}

#only subsets data with where the locus column of position_parsed matches the column name of gdata
#subsets every single cell in corresponding matching columns to return variant amino acids
#for(f in 1:length(position_parsed)){
# for(g in 1:length(gdata[which(colnames(gdata)%in%position_parsed[[f]][,1])])){
#  for(h in 1:nrow(gdata)){
#   variant_match[[g]][[h]][[f]]<-unique(subset(position_parsed[[f]][3], position_parsed[[f]][2]==gdata[which(colnames(gdata)%in%position_parsed[[f]][,1])][g][h,]))}
#}}




#removes elements with zero rows (i.e. no matches)
#for(f in 1:length(position_parsed)){
# for(g in 1:length(gdata[which(colnames(gdata)%in%position_parsed[[f]][,1])])){
#  for(h in 1:nrow(gdata)){
#   variant_match[[g]][[h]]<-variant_match[[g]][[h]][sapply(variant_match[[g]][[h]], nrow)>0]}}}
