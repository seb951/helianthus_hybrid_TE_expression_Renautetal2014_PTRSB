#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from = as.numeric(args[1])
to = as.numeric(args[2])


wd = "~/Documents/transposon"
setwd(wd) ### setup the working directory in the right format with the proper subdirectories.

################
### create file for each of the 7 comparisons ###
################
all = read.delim("reference/hybrid_species_new_2013-09-11", header = T, stringsAsFactors = F)
#snp_table = as.matrix(read.delim("results3/snp_table_all2_1to20000", header = T, sep = " ", stringsAsFactors = F))

#for(pp in from:to){

if(pp == "all") {snp_table = as.matrix(read.delim("results4/snp_table_all2_1_51k", header = T, sep = " ", stringsAsFactors = F)); snp_table[,2] = gsub(" ","",snp_table[,2]);
#change names
comp = snp_table ;
all = read.delim("reference/hybrid_species_new_2013-09-11", header = T, stringsAsFactors = F); colnames(comp)[4:ncol(comp)] = paste(colnames(comp)[4:ncol(comp)],all[,2],sep = "_");
snp_table = comp;


#create a specific dataset for each comparison? MAYBE NOT NECESSARY
ann_pet =  cbind(snp_table[,1:3],snp_table[,regexpr("_ann|_pet", colnames(snp_table))> 0]);
ann_des = cbind(snp_table[,1:3],snp_table[,regexpr("_ann|_des", colnames(snp_table))> 0]);
ann_ano = cbind(snp_table[,1:3],snp_table[,regexpr("_ann|_ano", colnames(snp_table))> 0]) ;
ann_par = cbind(snp_table[,1:3],snp_table[,regexpr("_ann|_par", colnames(snp_table))> 0]);
pet_des = cbind(snp_table[,1:3],snp_table[,regexpr("_pet|_des", colnames(snp_table))> 0]);
pet_ano = cbind(snp_table[,1:3],snp_table[,regexpr("_pet|_ano", colnames(snp_table))> 0]) ;
pet_par = cbind(snp_table[,1:3],snp_table[,regexpr("_pet|_par", colnames(snp_table))> 0]);

#
write.table(ann_pet, "results4/snp_table_ann_pet", row.names = F, col.names = T, quote = F);
write.table(ann_des, "results4/snp_table_ann_des", row.names = F, col.names = T, quote = F);
write.table(ann_ano, "results4/snp_table_ann_ano", row.names = F, col.names = T, quote = F);
write.table(ann_par, "results4/snp_table_ann_par", row.names = F, col.names = T, quote = F);
write.table(pet_des, "results4/snp_table_pet_des", row.names = F, col.names = T, quote = F);
write.table(pet_ano, "results4/snp_table_pet_ano", row.names = F, col.names = T, quote = F);
write.table(pet_par, "results4/snp_table_pet_par", row.names = F, col.names = T, quote = F)}

if(pp == 1) comp = as.matrix(read.delim("results4/snp_table_ann_pet", header = T, sep = " ", stringsAsFactors = F))
if(pp == 2) comp = as.matrix(read.delim("results4/snp_table_ann_des", header = T, sep = " ", stringsAsFactors = F))
if(pp == 3) comp = as.matrix(read.delim("results4/snp_table_ann_ano", header = T, sep = " ", stringsAsFactors = F))
if(pp == 4) comp = as.matrix(read.delim("results4/snp_table_ann_par", header = T, sep = " ", stringsAsFactors = F))
if(pp == 5) comp = as.matrix(read.delim("results4/snp_table_pet_des", header = T, sep = " ", stringsAsFactors = F))
if(pp == 6) comp = as.matrix(read.delim("results4/snp_table_pet_ano", header = T, sep = " ", stringsAsFactors = F))
if(pp == 7) comp = as.matrix(read.delim("results4/snp_table_pet_par", header = T, sep = " ", stringsAsFactors = F))


###
pairs = c("snp_table_ann_pet","snp_table_ann_des","snp_table_ann_ano","snp_table_ann_par","snp_table_pet_des","snp_table_pet_ano","snp_table_pet_par","snp_table_all2")

#comp = read.table(paste("results3/",pairs[pp],sep = ""), header = T, stringsAsFactors = F)
#comp = comp[c(1:2000),]

too_much_missing = c(1:nrow(comp))


	for(i in 1:nrow(comp))
		{
		n_ind = (ncol(comp)-3)
		too_much_missing[i] = length(c(1:n_ind)[grepl("XX",comp[i,4:ncol(comp)]) == T]) / n_ind
		if(i %% 100000 == 0) print(paste(i,Sys.time()))
		}

if(pp == 1) comp2 = comp[too_much_missing < 0.2,]
comp2 = comp[too_much_missing < 0.2,]
#if(pp == 2) comp2 = comp[too_much_missing < 0.3,]
#if(pp == 3) comp2 = comp[too_much_missing < 0.2,]



##############################
### trim based on 2pq (Ht) ###
##############################
comp3 =  c(1)

#counting the alleles.
	counter_ACGTX = matrix(0,nrow = nrow(comp2), ncol = 5)
	colnames(counter_ACGTX) = c("A","C","G","T","X")
	ht = c(1:nrow(comp2))
	n_ind = ncol(comp2)

	for(i in 1:nrow(comp2))
		{
		x = paste(comp2[i,4:n_ind], collapse = "")
		xx = strsplit(x, split = "")
		
		counter_ACGTX[i,1] = length(c(1:((n_ind-3)*2))[grepl("A", xx[[1]])])
		counter_ACGTX[i,2] = length(c(1:((n_ind-3)*2))[grepl("C", xx[[1]])])
		counter_ACGTX[i,3] = length(c(1:((n_ind-3)*2))[grepl("G", xx[[1]])])
		counter_ACGTX[i,4] = length(c(1:((n_ind-3)*2))[grepl("T", xx[[1]])])
		counter_ACGTX[i,5] = length(c(1:((n_ind-3)*2))[grepl("X", xx[[1]])])
		p = sort(counter_ACGTX[i,1:4])[4] 
		q = sort(counter_ACGTX[i,1:4])[3] 
		al = sum(counter_ACGTX[i,1:4])
		ht[i] = 2 * (p / al) * (q / al)
		}
	comp3 = comp2[ht > 0.2,]

#comp2 = comp = NULL

###################################
### trim based on Ho (paralogs) ###
###################################
comp4 = c(1)
ho  = rep(0,nrow(comp3))

for(i in 1:nrow(comp3))
	{
		a1 = substring(comp3[i,4:ncol(comp3)],1,1)
		a2 = substring(comp3[i,4:ncol(comp3)],2,2)		
		for(h in 1:length(a1))
		{
		if((a1[h] != a2[h]) & (a1[h] != "X") & (a2[h] != "X") & !is.na(a1[h]) & !is.na(a2[h])) ho[i] = (ho[i] + 1) #count the heterozygotes
		}
	}
	ho = ho / (ncol(comp3)-3) # observed heterozygosity
	comp4 = comp3[ho < 0.6,]

print(dim(comp4))

########################
### hierfstat format ###
########################
#takes about 1min per 1000 snps. 
library(hierfstat)

fstat_results = cbind(c(1:nrow(comp4)),0,0,0,0,0,0,0)

colnames(fstat_results) = c("snp","fst.ann.pet","fst.ann.des","fst.ann.ano","fst.ann.par","fst.pet.des","fst.pet.ano","fst.pet.par")

if(pp != "all") {

	for(i in 1:nrow(fstat_results))
	{
	hf = t( rbind( rep("empty",length( comp4[i,4:ncol(comp4)])),comp4[i,4:ncol(comp4)]))


	if(pp == all) { hf[regexpr("_ann",rownames(hf)) >0,1] = "1";  hf[regexpr("_pet",rownames(hf)) >0,1] = "2";  hf[regexpr("_des",rownames(hf)) >0,1] = "3";  hf[regexpr("_ano",rownames(hf)) >0,1] = "4";  hf[regexpr("_par",rownames(hf)) >0,1] = "5"}
	if(pp == 1) { hf[regexpr("_ann",rownames(hf)) >0,1] = "1";  hf[regexpr("_pet",rownames(hf)) >0,1] = "2"}
	if(pp == 2) { hf[regexpr("_ann",rownames(hf)) >0,1] = "1";  hf[regexpr("_des",rownames(hf)) >0,1] = "2"}
	if(pp == 3) { hf[regexpr("_ann",rownames(hf)) >0,1] = "1";  hf[regexpr("_ano",rownames(hf)) >0,1] = "2"}
	if(pp == 4) { hf[regexpr("_ann",rownames(hf)) >0,1] = "1";  hf[regexpr("_par",rownames(hf)) >0,1] = "2"}
	if(pp == 5) { hf[regexpr("_pet",rownames(hf)) >0,1] = "1";  hf[regexpr("_des",rownames(hf)) >0,1] = "2"}
	if(pp == 6) { hf[regexpr("_pet",rownames(hf)) >0,1] = "1";  hf[regexpr("_ano",rownames(hf)) >0,1] = "2"}
	if(pp == 7) { hf[regexpr("_pet",rownames(hf)) >0,1] = "1";  hf[regexpr("_par",rownames(hf)) >0,1] = "2"}
	
	hf  = gsub("U|R|Y|M|K|S|W|B|D|H|V|N","X",hf)
	hf = gsub("X","",hf)
	hf = gsub("0","",hf)
	hf = gsub("A",1, hf)
	hf = gsub("C",2, hf)
	hf = gsub("G",3, hf)
	hf = gsub("T",4, hf)

	colnames(hf) = c("pop","snp1")
		s = rbind(c(1,2),c(1,3),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5)) #7 species pairs for which you need to calculate fst
	
		for(ii in 1:7)
			{
				hf1 = hf[hf[,1] == s[ii,1], ]; hf1 = rbind(hf1,hf[hf[,1] == s[ii,2], ])
				if(length(hf1[hf1[,2] == "",2]) /length(hf1[,2]) < 0.5) x = varcomp(data.frame(as.integer(hf1[,1]), as.integer(hf1[,2])))
				if(!is.na(x$F[1,1]))	fstat_results[i,ii+1] = signif(x$F[1,1],4) #Fst
			}
		if(i %% 1000 == 0) print(paste(i,"of",nrow(fstat_results),Sys.time()))
	}
	}

comp4 = cbind(comp4, fstat_results[,1:2])
####################
if(pp == "all") write.table(comp4,"results4/comp4_allspecies", col.names = T, row.names = F, quote = F)
if(pp == 1) write.table(comp4,"results4/comp4_annpet", col.names = T, row.names = F, quote = F)
if(pp == 2) write.table(comp4,"results4/comp4_anndes", col.names = T, row.names = F, quote = F)
if(pp == 3) write.table(comp4,"results4/comp4_annano", col.names = T, row.names = F, quote = F)
if(pp == 4) write.table(comp4,"results4/comp4_annpar", col.names = T, row.names = F, quote = F)
if(pp == 5) write.table(comp4,"results4/comp4_petdes", col.names = T, row.names = F, quote = F)
if(pp == 6) write.table(comp4,"results4/comp4_petano", col.names = T, row.names = F, quote = F)
if(pp == 7) write.table(comp4,"results4/comp4_petpar", col.names = T, row.names = F, quote = F)

if(pp == "all_dont to it") {
comp4 = read.table("results4/comp4_allspecies", header = T, stringsAsFactors = F);
comp4_ptgs = read.table("results4/comp4_ptgs", header = T, stringsAsFactors = F);
comp4_transpos = read.table("results4/comp4_transpos", header = T, stringsAsFactors = F);

### ###	###	###
###	plotting ###	###
###	###	###	###
###	are fst higher for PTGS genes? ### ANSWER : NO, THEY ARE THE SAME!

par(mfrow = c(3,3));

for(ii in 1:7)
{
breaks = seq(-1.1,1.1,by = 0.01)
x = hist(comp4[comp4[,48+ii] != 0,48+ii], breaks = breaks,plot = F)
plot(x$mids[90:220],x$density[90:220], type = "h",lwd = 5,xlab = colnames(comp4)[48+ii], col = "#00000050")

breaks = seq(-1.1,1.1,by = 0.01)
xx = hist(comp4_ptgs[comp4_ptgs[,48+ii] != 0,48+ii], breaks = breaks,plot = F)

points(xx$mids[90:220],xx$density[90:220], type = "h",lwd = 5,xlab = colnames(comp4_ptgs)[48+ii], col = "#00009950")
}
}
#dev.print(device=svg, "~/Desktop/distri_null.svg", onefile=FALSE)
}


