#!/usr/bin/Rscript --verbose

###Step 0.1 set up working directory ###
setwd("/home/seb/SITES/transposon/")
######################
###get the length of the genes###
###################### 
ref = read.delim("/home/seb/Documents/transposon/reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F, stringsAsFactors = F)
x = c(1:nrow(ref))[gregexpr(">",ref[,1]) > 0]
len = matrix(0, nrow = length(x), ncol = 7); len[,1] =  ref[x,1]; x = c(x,nrow(ref))

colnames(len)  = c("name","len", "pi_ann", "pi_ano","pi_des","pi_par","pi_pet")

for(i in 1:nrow(len)) {len[i,2] = sum(nchar(ref[(x[i]+1):(x[i+1]-1),1]))}
len[,1] = gsub(">","", len[,1],fixed = T)
###

comp4 = read.delim("/home/seb/Documents/transposon/results4/comp4_allspecies", header = T,stringsAsFactors = F, sep =  " ")
all = read.delim("/home/seb/Documents/transposon/reference/hybrid_species_new_2013-09-11", header = T, stringsAsFactors = F); colnames(comp4)[4:ncol(comp4)] = paste(colnames(comp4)[4:ncol(comp4)],all[,2],sep = "_")


###
aaa = Sys.time()
un_comp4 = unique(comp4[,1]) # unique set of genes, this is just faster than looking through the whole matrix everytime.
	for(i in 1:nrow(len)) 
#	for(i in 1:1000)
	{
	if(length(un_comp4[un_comp4 == len[i,1]]) > 0)  {
	
	x1 = comp4[comp4[,1] == len[i,1],4:57];
	
	if(is.na(x1[1,1]))	{x1 = matrix("AA", nrow = 2, ncol = 54)}; # just make sure the data is in a matrix format, even if nothing is there. 
	if(nrow(x1) == 1)	x1 = rbind(x1,"AA"); # if only one row, add another one to treat as matrix.

	x1 =  gsub("0","XX",as.matrix(x1));	
	x1 =  gsub("X","N",as.matrix(x1));	
	
	x1_one = x1; #x1[,regexpr(sps[o], colnames(x1))> 0]
	one_presites = matrix(0,nrow = (ncol(x1_one)*2), ncol = 2);
	seq = NULL; for(r in 1:54) {seq = c(seq,rep(r,2))}; one_presites[,1] = substring(all[seq,2],1,3); #
	for(n in 1:ncol(x1_one))
	{one_presites[(n*2)-1,2] = paste(substring(x1_one[,n],1,1),collapse = "");one_presites[(n*2),2] = paste(substring(x1_one[,n],2,2),collapse = "")};
	one_presites[,2] =  paste(one_presites[,1],paste(one_presites[,2],paste(rep("A",as.numeric(len[i,2])-nrow(x1_one)),collapse = ""),sep = ""),sep = "       ");
	one_presites_ordered = one_presites[order(one_presites[,1]),];
	
	#write tex_deb file 
	line1 = "SITES sample input" ;
	line2 = paste(ncol(x1_one)*2 -16, nchar(one_presites[1,2])-9); ###16 is to remove the FUCKING f1s...
	line3 = paste(1,1);
	line4 = paste(1,nchar(one_presites[1,2])-1);
 	line5	= 5;
	line6 = paste(rle(one_presites_ordered[,1])$values[1],rle(one_presites_ordered[,1])$lengths[1]);
	line7 = paste(rle(one_presites_ordered[,1])$values[2],rle(one_presites_ordered[,1])$lengths[2]);
	line8 = paste(rle(one_presites_ordered[,1])$values[3],rle(one_presites_ordered[,1])$lengths[3]);
	line9 = paste(rle(one_presites_ordered[,1])$values[5],rle(one_presites_ordered[,1])$lengths[5]); ###skip one here to remove the f1s...
	line10 = paste(rle(one_presites_ordered[,1])$values[6],rle(one_presites_ordered[,1])$lengths[6]);
	
	sites = c(line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,one_presites[one_presites[,1] != "f1",2]);

	write.table(sites, "sites_sample", row.names = F, col.names = F, quote = F);

	system("./sites -isites_sample -rsites_results -mmy_data -sa -asp -ogc -cf");

###read output###
#system("grep 'Fst' sites_results.SIT -A 5 | tail -1 >fst") 
	system("grep 'Theta' sites_results.SIT -A 6 | tail -5 >pi") ;

#fst = as.numeric(read.table("fst")[3])
	pi = read.table("pi")[8];

###update the table 
	len[i,3:7] = pi[,1]}
	
	
	
		if(i %% 1000 == 0) print(paste(i,Sys.time()))
	}
Sys.time() - aaa

#write.table(len,"/home/seb/Documents/helianthus_analysis/92samples_new_alignments_sept11/results/fst_pergene", row.names = F, col.names =T)
write.table(len,"SITES_pi", row.names = F, col.names =T)




