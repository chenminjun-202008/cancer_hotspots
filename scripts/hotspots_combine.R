args <- commandArgs(T)

if (length(args) != 2){

		print("Rscript.R <script.R> <hostspot.xls> <out.txt>")
		q()
}

#--------------
library(tidyverse)
library(readxl)


# 读取信息表并将2个sheet合并

data1 <- read_excel(args[1], sheet = "SNV-hotspots") %>% unique() 
data2 <- read_excel(args[1], sheet = "INDEL-hotspots") %>% unique()


data1 <- data1 %>% 
	extract(Reference_Amino_Acid,remove=FALSE, regex="^(.*):(.*)$",into=c("ref_AA","pos1"))%>%
	extract(Variant_Amino_Acid,remove=FALSE, regex="^(.*):(.*)$",into=c("alt_AA","positive_num"))%>%
	mutate(Amino_Acid_Change=paste0(ref_AA,Amino_Acid_Position,alt_AA)) %>%
	mutate(Amino_Acid_Change=ifelse(
			grepl("splice",Amino_Acid_Position),
			gsub("_splice","_sice",Amino_Acid_Position),
			Amino_Acid_Change)) %>%
	mutate(Amino_Acid_Change=paste0(Hugo_Symbol,":",Amino_Acid_Change)) %>%
	select(Amino_Acid_Change,positive_num,everything()) %>% 
	select(-c(ref_AA,pos1,alt_AA))

data2 <- data2 %>% 
	extract(Variant_Amino_Acid,remove=FALSE, regex="^(.*):(.*)$",into=c("Amino_Acid_Change","positive_num")) %>%
	mutate(Amino_Acid_Change=paste0(Hugo_Symbol,":",Amino_Acid_Change)) 


#-----------------------
header <- unique(c(colnames(data1),colnames(data2)))

for (x in header[!header %in% names(data1)]) { data1[,x] <- "" }
for (x in header[!header %in% names(data2)]) { data2[,x] <- "" }

#colnames(data1)
#colnames(data2)
merge <- rbind(data1,data2) %>%
	write_tsv(args[2])
	



