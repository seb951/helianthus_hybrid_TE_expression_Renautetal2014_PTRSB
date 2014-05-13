#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from_cmd = as.numeric(args[1]) # Specify which sequences in "list_ind" file you want to align, directlty from the shell. Alternatively, you can do this from the "alignments" function itself.
to_cmd = as.numeric(args[2])


#setwd("~/Documents/transposon") ### setup the working directory in the right format
setwd("~/transposon")

trimo = function(list_ind= "reference/hybrid_species", from = 1, to = 1) {

individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals)) {individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])]} # proper format

for(i in from : to)
{
###trimommatic #####
phred = ifelse(individuals[i,4] == "sanger","-phred33","-phred64")
seq_1_in = gsub(".fq","_1.fq",individuals[i,1])
seq_2_in = gsub(".fq","_2.fq",individuals[i,1])
seq_1_out_p = gsub(".fq","_trim1.fq",individuals[i,5]);seq_1_out_p = sub("^","data/",seq_1_out_p)
seq_2_out_p = gsub(".fq","_trim2.fq",individuals[i,5]);seq_2_out_p = sub("^","data/",seq_2_out_p)
seq_1_out_up = gsub(".fq","_unpaired1.fq",individuals[i,5]);seq_1_out_up = sub("^","data/",seq_1_out_up)
seq_2_out_up = gsub(".fq","_unpaired2.fq",individuals[i,5]);seq_2_out_up = sub("^","data/",seq_2_out_up)
log = gsub(".fq","_trimo.log",individuals[i,5]);log = sub("^","data/",log)

trim = paste("java -classpath ~/Trimmomatic-0.25/trimmomatic-0.25.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 12",phred,"-trimlog",log,seq_1_in,seq_2_in,seq_1_out_p,seq_1_out_up,seq_2_out_p,seq_2_out_up,"ILLUMINACLIP:reference/IlluminaContaminants.fa:2:40:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:10:15 MINLEN:36",sep = " ")

system(trim)
system(paste("rm", seq_1_out_up, seq_2_out_up, log))
}
}


###################
### run TRIMMOMATIC ###
###########e########
trimo(list_ind= "reference/hybrid_species", from = from_cmd, to = to_cmd)


