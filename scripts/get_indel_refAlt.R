args <- commandArgs(T)

if(length(args) != 4){
	print("Rscript script.R <hotsplot.txt> <input.txt> <input.trasvar> <outfile>")
	q()

}

#================
library(tidyverse)
out <- args[4]

transvar <- read_tsv(args[3]) %>% select(input,CHROM,POS,REF,ALT) %>% unique %>%
	mutate(CHROM=gsub("chr","",CHROM))
transvar[] <- lapply(transvar,as.character)

input <- read_tsv(args[2],col_names=c("input","Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2")) %>% unique
input[] <- lapply(input,as.character)

result <- input %>%
	left_join(transvar,by=c("input")) %>%
	select(Chromosome:Tumor_Seq_Allele2,CHROM:ALT) %>% unique %>%
	write_tsv("tmp.txt")

data <- read_tsv(args[1])
data[] <- lapply(data,as.character)

#-------
snv <- data %>% filter(Reference_Allele!="-" & Tumor_Seq_Allele2!="-")

indel <- data %>% filter(!(Reference_Allele!="-" & Tumor_Seq_Allele2!="-")) %>%
	left_join(result,by=c("Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2")) %>%
	unique
check <- indel %>% filter(is.na(CHROM)) %>% write_tsv("need_check.txt")
message("Need check",nrow(check))

indel2 <- indel %>% mutate(Chromosome=CHROM,Start_Position=POS,Reference_Allele=REF,Tumor_Seq_Allele2=ALT) %>%
	select(-c(CHROM:ALT))

#--
rbind(snv,indel2) %>% write_tsv(out) 



