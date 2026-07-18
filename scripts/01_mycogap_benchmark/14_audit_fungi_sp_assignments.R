# Purpose: Audit ambiguous kingdom-only fungal assignments produced by fungi-only reference databases and classifiers.


library(dada2)
library(data.table)
library(tidyverse)




thread <- 16
minBoot <- 50

# use the HMP, WeP1, and CHGM data as examples
seqtab_path <- "/data/wangxinyu/Fungi_Atlas/Data/Public/Annotation/PE/PRJNA356769b_ITS2/1_dada2/result/PRJNA356769b_ITS2_seqtab.csv"
seqtab_path <- "/data/wangxinyu/Fungi_Atlas/Data/WeGut/WeP1_ITS2/1_dada2/result/WeP1_ITS2_seqtab.csv"
seqtab_path <- "/data/wangxinyu/Fungi_Atlas/Data/Public/Annotation/PE/PRJCA010668_ITS1/1_dada2/result/PRJCA010668_ITS1_seqtab.csv"

# use the fungi-only ref dataset (https://doi.org/10.15156/BIO/3301229)
ref <- "/data/wangxinyu/Ref/ITS_anno/Unite/sh_general_release_19.02.2025/sh_general_release_dynamic_19.02.2025.fasta"

# start taxonomy assignment
seqtab.nochim <- fread(seqtab_path, data.table = FALSE)
seqtab.nochim <- seqtab.nochim %>% column_to_rownames("V1")
seqtab.nochim[1:3, 1:3]
seqs <- colnames(seqtab.nochim)
length(seqs)
head(seqs)

batch_size <- thread * 100
batch_size
num_batches <- ceiling(length(seqs) / batch_size)
num_batches

taxa_list <- list()
boot_list <- list()

for (i in 1:num_batches) {
  up <- min(i * batch_size, length(seqs))
  down <- (i - 1) * batch_size + 1
  seqs_batch <- seqs[down:up] # Select sequences for this batch
  taxa_batch <- assignTaxonomy(seqs_batch, ref, minBoot = minBoot, outputBootstraps = TRUE, multithread = thread, tryRC = TRUE, verbose = TRUE)

  taxa_list[[i]] <- taxa_batch[["tax"]]
  boot_list[[i]] <- taxa_batch[["boot"]]

  rm(seqs_batch, taxa_batch) # Free up memory
  cat(paste0("\n", "-- batch", i, " / ", num_batches,": seq", down, "-", up ," was down! --", "\n"))
}

taxa <- do.call(rbind, taxa_list)
boot <- do.call(rbind, boot_list)


# check whether there is the "fungi_sp"
taxa <- taxa %>% as.data.frame()
boot <- boot %>% as.data.frame()

unique(taxa$Kingdom)
unique(boot$Kingdom)
# HMP, WeP1, and CHGM
# [1] "k__Fungi"
# [1] 100

unique(taxa$Phylum)
# HMP
#[1] "p__Ascomycota"    "p__Basidiomycota" NA                 "p__Mucoromycota"

# WeP1
# [1] "p__Ascomycota"               NA
# [3] "p__Basidiomycota"            "p__Mucoromycota"
# [5] "p__Fungi_phy_Incertae_sedis" "p__Rozellomycota"
# [7] "p__Chytridiomycota"          "p__Mortierellomycota"
# [9] "p__Neocallimastigomycota"

# CHGM
# [1] "p__Ascomycota"               "p__Basidiomycota"
# [3] NA                            "p__Mortierellomycota"
# [5] "p__Neocallimastigomycota"    "p__Chytridiomycota"
# [7] "p__Glomeromycota"            "p__Fungi_phy_Incertae_sedis"
# [9] "p__Mucoromycota"             "p__Rozellomycota"
# [11] "p__Olpidiomycota"            "p__Kickxellomycota"
# [13] "p__Entorrhizomycota"         "p__Blastocladiomycota"
# [15] "p__Zoopagomycota"


sum(is.na(taxa$Phylum)) / nrow(taxa)

# HMP: 0.2601385
# WeP1: 0.6134499
# CHGM: 0.121026

# Write output
path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/fungi_sp"
write.csv(taxa, file = file.path(path, "dada2_taxa_hmp.csv"))
write.csv(boot, file = file.path(path, "dada2_boot_hmp.csv"))
write.csv(taxa, file = file.path(path, "dada2_taxa_wep1.csv"))
write.csv(boot, file = file.path(path, "dada2_boot_wep1.csv"))
write.csv(taxa, file = file.path(path, "dada2_taxa_chgm.csv"))
write.csv(boot, file = file.path(path, "dada2_boot_chgm.csv"))





q2_tax <- read.delim("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/fungi_sp/qiime2/result/taxtab_sklearn/taxonomy.tsv",
  header = TRUE,
  sep = "\t",
  quote = "")

head(q2_tax)


(q2_tax %>% filter(Taxon == "k__Fungi") %>% nrow()) / (q2_tax %>% nrow()) # 0.3307626
