# Purpose: Summarize the proportions of putatively novel ITS/16S features and their assigned reads.

library(tidyverse)


# data in
blast_ITS <- read.delim("/data/wangxinyu/Fungi_Atlas/Analysis/Main/blast/fungi/refseq_filterdp_derep_unite_final.tsv", header = TRUE)

# clean
blast_ITS <- blast_ITS %>% select (qseqid, pident)
blast_ITS$category <- "ITS"
head(blast_ITS)
nrow(blast_ITS) # 55902


per1_ITS <- round(mean(blast_ITS$pident < 97) * 100, 1)
per1_ITS # 31.4



ASV_counts0_ITS <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Main/data/ASV/refseq_filterdp_derep_counts.csv", row.names = 1)
head(ASV_counts0_ITS)
nrow(ASV_counts0_ITS) # 57197


blast_ITS$qseqid <- gsub(";size=\\d+$", "", blast_ITS$qseqid)
head(blast_ITS)
nrow(blast_ITS)

ASV_counts_ITS <- merge(blast_ITS, ASV_counts0_ITS, by.x = "qseqid", by.y = "row.names", all.x = TRUE)
head(ASV_counts_ITS)



count2_ITS <- ASV_counts_ITS %>%
  filter(pident < 97) %>%
  summarise(sum_counts = sum(counts, na.rm = TRUE)) %>%
  pull(sum_counts)
count2_ITS

per2_ITS <- round(count2_ITS / sum(ASV_counts_ITS$counts) * 100, 1)
per2_ITS # 3.9


# 16S
blast_16S <- read.delim("/data/wangxinyu/Fungi_Atlas/Analysis/Main/blast/cell2025/analysis/cell2025_asv_filter2_0.05_silva_final.tsv", header = TRUE)

# clean
blast_16S <- blast_16S %>% select (qseqid, pident)
blast_16S$category <- "16S"
head(blast_16S)
nrow(blast_16S)


per1_16S <- round(mean(blast_16S$pident < 97) * 100, 1)
per1_16S # 17.4



ASV_counts0_16S <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Main/blast/cell2025/analysis/cell2025_asv_filter_counts.csv", row.names = 1)
head(ASV_counts0_16S)
nrow(ASV_counts0_16S)


ASV_counts_16S <- merge(blast_16S, ASV_counts0_16S, by.x = "qseqid", by.y = "ASV", all.x = TRUE)
head(ASV_counts_16S)
nrow(ASV_counts_16S)



count2_16S <- ASV_counts_16S %>%
  filter(pident < 97) %>%
  summarise(sum_counts = sum(counts, na.rm = TRUE)) %>%
  pull(sum_counts)

count2_16S

per2_16S <- round(count2_16S / sum(ASV_counts_16S$counts) * 100, 1)
per2_16S # 0.7


# overall
per1_ITS # 31.4
per2_ITS # 3.9
per1_16S # 17.4
per2_16S # 0.7
