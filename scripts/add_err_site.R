args <- commandArgs(T)

if(length(args) != 3){
	print("Rscript SCRIPT.r <hotspots.txt> <transvar.txt> <final.txt>")
	q()
}

#=================
library(tidyverse)

hot <- read_tsv(args[1])
transvar <- read_tsv(args[2])
out <- args[3]

#-------#chr2:g.26101089C>G/c.3G>C/p.M1I
transvar <- transvar %>% select(input,transcript,"coordinates(gDNA/cDNA/protein)") %>%
	rename(mutations="coordinates(gDNA/cDNA/protein)") %>%
	filter(transcript!=".",input!="input") %>%
	mutate(input=str_split(input,"\\|",simplify=T)[,2]) %>%
	extract(mutations,remove=FALSE, regex="^(.*)/(.*)/(.*)$",into=c("gDNA","cDNA","protein")) %>%
	extract(gDNA,remove=FALSE, regex="^(.*):g.(.*)(.)>(.)$",into=c("chrom","pos","refs","alts")) %>%
	select(input,mutations,chrom,pos,refs,alts,cDNA) %>%
	unique %>% write_tsv("tmp.txt")

colnames(hot)
colnames(transvar)
hot <- hot %>% left_join(transvar,by=c("Amino_Acid_Change"="input")) %>%
	mutate(Chromosome=gsub("chr","",chrom),
				 Start_Position=pos,
				 End_Position=pos,
				 Reference_Allele=refs,
				 Tumor_Seq_Allele2=alts,
				 HGVSc=cDNA) %>%
	select(-c(mutations:alts,cDNA)) %>%
	unique %>%
	write_tsv(out)
