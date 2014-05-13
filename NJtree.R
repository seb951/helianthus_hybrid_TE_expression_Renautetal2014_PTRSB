#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
species = args[1] # can specify which species you want to align, directlty from the shell. Alternatively, you can do this from the "snp" function itself.
cpu = as.numeric(args[2])
################
###requirements #####
################

### YOU NEED SAMTOOLS and BCFTOOLS. You should also have vcfutils.pl in your "Rcode/" directory ###
#you should run snp() species per species, but not simultaneously. Otherwise there will be some problems with overwriting files. Uncomment the 
#you should run snp_stats() once you have ran snp for all 3 species.

###########################
### defining the SNP calling function ###
###########################
splitstree = function(list_ind= "reference/hybrid_species",user = "renaut", species = "all",genes = 1000,cpu = 16,wd = "~/transposon") {

setwd(wd) ### setup the working directory in the right format with the proper subdirectories.

###create "split" object for mpileup if it doesn't exist###
system("grep '>' reference/HA412_trinity_noAltSplice_400bpmin.fa | awk '{print $1}' >reference/temp")
temp = read.table("reference/temp", stringsAsFactors = F, header = F);temp[,1] = cbind(gsub(">","",temp[,1]),"A")

if(file.exists("reference/split") == F){
for(i in 1: nrow(temp)) {
cmd = paste("awk 'NR == ", i*2,"' reference/HA412_trinity_noAltSplice_400bpmin.fa | wc -m >reference/temp_wc", sep = ""); system(cmd)
temp[i,2] = paste("1-",read.table("reference/temp_wc",stringsAsFactors = F, header = F)-1,sep = "")
if(i %% 10000 == 0) print(paste(i,"of 51k, Time is:",Sys.time()))
}
split = as.data.frame(paste(temp[,1],temp[,2],sep = ":"));colnames(split) = "V1";rm(temp) #split file for mpileup
} else split = read.table("reference/split", stringsAsFactors = F)

###load individuals' names
individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals))
{individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])];
if(length(grep("fq",individuals[p,1])) == 1) {individuals[p,5] = paste(individuals[p,5],".sorted.bam",sep = ""); individuals[p,1] = paste(individuals[p,1],".sorted.bam",sep = "")}
if(length(grep("fq",individuals[p,1])) != 1) {individuals[p,5] = paste(individuals[p,5],".fq.sorted.bam",sep = ""); individuals[p,1] = paste(individuals[p,1],".fq.sorted.bam",sep = "")}
} # proper format to further process with idxstats

#all_final = paste("alignments/", individuals[,5], collapse = " ", sep = "")
all_final = paste(individuals[,1], collapse = " ", sep = "")

############
### mpileup ###
############
#To run mpileup gene per gene, this way you can parallelize.#
for(j in 1:genes)#2000 genes only
{
system(paste("ps -u",user,"| grep 'samtools' | wc -l >prog"))
mpileup = paste("samtools mpileup -C50 -I -ugf reference/HA412_trinity_noAltSplice_400bpmin.fa ", all_final, " -r ",split[j,1]," | bcftools view -cvg - > mpileup/variants",j,".raw.vcf", sep = "")
Sys.sleep(ifelse(read.table("prog") <= cpu,0, (read.table("prog") - cpu)^4    ))
mpileup_exe = paste("mpileup/cmd/mpileup_",j,sep = "")
write.table(mpileup,mpileup_exe, row.names = F, col.names = F, quote = F)
system(paste("chmod +x",mpileup_exe))
system(paste("nohup ./", mpileup_exe, ">log&",sep = ""))
if(j %% 100 == 0) print(paste("Calling SNPs of gene number",j,"(out of  ",genes, ") The time is:",Sys.time()))
}

#cat the variants_j.raw.vcf files. 
Sys.sleep(60) #tidy up the processes running, once you are done
system("echo -n '' > mpileup/0_51k_variants.raw.vcf") #final file
system("cat mpileup/variants*  >>mpileup/0_51k_variants.raw.vcf") # cat verything
system("grep '#' -v  mpileup/0_51k_variants.raw.vcf >mpileup/0_51k_variants.clean.vcf") #tidy up
system("rm mpileup/variants*")
system("rm mpileup/cmd/*")


###Pre-cleaning step###
	system("wc -l mpileup/0_51k_variants.clean.vcf >mpileup/wordcount1")
	wordcount = read.table("mpileup/wordcount1")
	system("awk 'NR == 27' mpileup/0_51k_variants.raw.vcf  >mpileup/header")
	header = read.delim("mpileup/header", header = F, stringsAsFactors = F)
	header = gsub("/home/seb/Documents/repeatability/","",header, fixed = T)
	header = gsub("/home/seb/Documents/repeatability/alignments/","",header, fixed = T)
	header = gsub("/media/seb_1TB/hybrid_species_MER_GRNnote_transposon/alignments/","",header, fixed = T)
	header = gsub("alignments/","",header[10:length(header)], fixed = T)
	header = gsub(".fq.sorted.bam","",header, fixed = T)
	header = gsub(".sorted.bam","",header, fixed = T)
	header = gsub("_R1.fq.bz2","",header, fixed = T)
	header = gsub("/SciBorg/array0/renaut/speciation_islands_individuals/","", header, fixed = T)
	write.table(t(as.matrix(c("reference", "position","reference_annuus_allele",header))),"results4/snp_table_all",row.names = F, col.names = F, quote = F, sep = " ")

	temp = NULL
	mis = NULL
	con = file("mpileup/0_51k_variants.clean.vcf")
	open(con)
	for(i in 1:as.numeric(wordcount[1]))
	{
		x = readLines(con,1) #you will now read the mega big 0-51k_variants.clean.vcf file line by line in R. 
		xx = strsplit(x, split = "\t")[[1]]
		xxx = xx[(length(xx)-(length(header) - 1)):length(xx)]
		
		for(q in 1:length(xxx))	# this is to replace low quality calls (maximum Phred-scaled genotype likelihoods below 20). ie. essentially, you need at least 2 reads to call a SNP!
		{
		temp = strsplit(xxx[q], split = ":|,")[[1]]
		if(max(as.numeric(temp[length(temp)])) < 20) xxx[q] = "XX" else {xxx[q] = temp[1]; xxx[q] = gsub("/","",xxx[q],fixed = T); xxx[q] = gsub("0","R",xxx[q],fixed = T); xxx[q] = gsub("1","A",xxx[q],fixed = T)} #new way of doing things...
		}
		xxx = gsub("R",xx[4],xxx)
		xxx = gsub("A",xx[5],xxx)
		xxx = c(xx[c(1,2,4)],xxx)
		if(nchar(xx[5]) > 1) xxx[4:length(xxx)] = rep("XX",length(xxx)-3)
		
		if((length(xxx[xxx == "XX"]) / length(xxx[4:length(xxx)])) < 0.4)	cat(t(as.matrix(xxx)),file = "results4/snp_table_all",append = T, fill = F, "\n") else mis = c(mis, i) #only cat loci with less than 40% missing data. 
	#	if(regexpr(xxx[3],paste(xxx[4:21],collapse = "")) < 0) temp = rbind(temp,xxx)
		if(i %% 1000 == 0) print(paste(i,"of", wordcount[1], Sys.time()))
		}
close(con)

system("cat results4/snp_table_all | sed 's/[ \t]*$//' >results4/snp_table_all2")
#system("rm mpileup/0_51* mpileup/header mpileup/wordcount1 log prog ") ### rm undesirable outputs
}

#########################
### RUN SNP CALLING FUNCTION ###
#########################
splitstree(list_ind= "reference/hybrid_species_new_2013-09-11",user = "seb",species = "all",cpu = 7, genes = 1000,wd = "~/Documents/transposon")


###############################
#### NEIGHBOUR JOINING TREE ###
###############################
snp_table = as.matrix(read.delim("results4/snp_table_all2", header = T, sep = " ", stringsAsFactors = F))
snp_table[,2] = gsub(" ","",snp_table[,2])
#individuals = read.delim("reference/all_species_nov2012_cleaned.txt", header = T, stringsAsFactors = F) #individuals of interest (454 individuals removed, weird individuals based on NatCom paper also removed#

################
### TRIMMING ###
################
all_comparisons = snp_table
###kick out SNP which have more than 20% XX data missing
too_much_missing = c(1:nrow(all_comparisons))

	for(i in 1:nrow(all_comparisons))
		{
		n_ind = ncol(all_comparisons)-3
		too_much_missing[i] = length(c(1:n_ind)[grepl("XX",all_comparisons[i,4:ncol(all_comparisons)]) == T]) / n_ind
		}

all_comparisons_2 = all_comparisons[too_much_missing < 0.1,]

##############################
### trim based on 2pq (Ht) ###
##############################
all_comparisons_3 = 1

#counting the alleles.
	counter_ACGTX = matrix(0,nrow = nrow(all_comparisons_2), ncol = 5)
	colnames(counter_ACGTX) = c("A","C","G","T","X")
	ht = c(1:nrow(all_comparisons_2))
	n_ind = ncol(all_comparisons_2)
	
	for(i in 1:nrow(all_comparisons_2))
	#for(i in 1:100)
		{
		x = paste(all_comparisons_2[i,4:n_ind], collapse = "")
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
	all_comparisons_3 = all_comparisons_2[ht > 0.2,]


###################################
### trim based on Ho (paralogs) ###
###################################

all_comparisons_4 = list(1)

ho  = rep(0,nrow(all_comparisons_3))

for(i in 1:nrow(all_comparisons_3))
	{
		a1 = substring(all_comparisons_3[i,4:ncol(all_comparisons_3)],1,1)
		a2 = substring(all_comparisons_3[i,4:ncol(all_comparisons_3)],2,2)		
		for(h in 1:length(a1))
		{
		if((a1[h] != a2[h]) & (a1[h] != "X") & (a2[h] != "X") & !is.na(a1[h]) & !is.na(a2[h])) ho[i] = (ho[i] + 1) #count the heterozygotes
		}
	}
	ho = ho / (ncol(all_comparisons_3)-3) # observed heterozygosity
	all_comparisons_4 = all_comparisons_3[ho < 0.6,]

###################################
### create a single consensus of both alleles. ###
###################################	
col = ncol(all_comparisons_4)
all_comparisons_4_nj_preformat = all_comparisons_4

for(i in 1:nrow(all_comparisons_4))
	{
	temp1 = all_comparisons_4[i,4:col]
	temp2 = temp1
	for(j in 1:length(temp2))
		{
		if(substring(temp1[j],1,1) == substring(temp1[j],2,2)) temp2[j] = substring(temp1[j],1,1)
		if(substring(temp1[j],1,1) == "X" & substring(temp1[j],2,2) != "X" ) temp2[j] = substring(temp1[j],2,2)
		if(substring(temp1[j],2,2) == "X" & substring(temp1[j],1,1) != "X" ) temp2[j] = substring(temp1[j],1,1)
	
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "C") | (substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "A"))) temp2[j] = "M"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "G") | (substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "A"))) temp2[j] = "T"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "A"))) temp2[j] = "W"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "G") | (substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "C"))) temp2[j] = "S"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "C"))) temp2[j] = "Y"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "G"))) temp2[j] = "K"
		
		if(temp2[j] == "X") temp2[j] = "N"
		}
	if((i %% 500) == 0) print(paste(i, "of",nrow(all_comparisons_4),"the time is now:",Sys.time()))
	all_comparisons_4_nj_preformat[i,4:col] = temp2
	}
#all_comparisons_4_nj_preformat = all_comparisons_4_nj_preformat[1:5000,]
write.table(all_comparisons_4_nj_preformat,"results4/all_comparisons_4_nj_preformat", row.names = F, col.names = T, quote = T)	

##############
####NEXUS FILE FOR SPLITSTREE###
##############
all_comparisons_4_nj_preformat = as.matrix(read.delim("results4/all_comparisons_4_nj_preformat", header  = T, sep = " "))
individuals = read.delim("reference/hybrid_species_new_2013-09-11", header = T, stringsAsFactors = F)

colnames(all_comparisons_4_nj_preformat)[4:57] = individuals[,5]

sequences = c(1:(ncol(all_comparisons_4_nj_preformat) - 3))


for(j in 4:ncol(all_comparisons_4_nj_preformat))
{
temp = colnames(all_comparisons_4_nj_preformat)[j]
temp = paste(temp,paste(all_comparisons_4_nj_preformat[,j], collapse = ""), sep = "     ") #create one consensus concatenated sequence per individual. 
sequences[j-3] = temp
}

nrow(all_comparisons_4_nj_preformat)
ncol(all_comparisons_4_nj_preformat)



all_helianthus_speciation_island = c("#NEXUS","BEGIN taxa;",paste("DIMENSIONS ntax=",ncol(all_comparisons_4_nj_preformat)-3,";", sep = ""),"TAXLABELS",colnames(all_comparisons_4_nj_preformat)[4:ncol(all_comparisons_4_nj_preformat)],";","END;","BEGIN characters;",paste("DIMENSIONS nchar=",nrow(all_comparisons_4_nj_preformat),";",sep = ""),"FORMAT","datatype=DNA","missing=N","gap=-","symbols=\"A C G T\"","labels","interleave",";","MATRIX",sequences, ";","END;") #CRAZY splitstree formating

write.table(all_helianthus_speciation_island,"results4/all_helianthus_hybrid_species_20k.nex", sep = "", quote = F, row.names = F,col.names = F)





############
###SANDBOX### 
############




