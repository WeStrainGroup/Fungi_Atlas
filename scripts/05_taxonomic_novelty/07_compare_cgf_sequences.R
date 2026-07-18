# Purpose: Extract CGF ITS1/ITS2 sequences, integrate atlas matches, and quantify culture-supported novelty.
# Panels: Figure S5G-H.

library(tidyverse)
library(Biostrings)
library(ggridges)



meta <- read.delim("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/raw/gut_fungi_meta.tsv")
head(meta)

gcf_id <- meta %>% filter(tag == "CGF_Catalog") %>% pull(genome_id)
head(gcf_id)
length(gcf_id) # 708

# input all seq
seq_ITS1 <- readDNAStringSet("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/raw/obipcr_ITS1_ITSx.fasta")
seq_ITS2 <- readDNAStringSet("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/raw/obipcr_ITS2_ITSx.fasta")
seq_ITS1
seq_ITS2


seq_ITS1_filter <- seq_ITS1[sub("__.*", "", names(seq_ITS1)) %in% gcf_id]
seq_ITS2_filter <- seq_ITS2[sub("__.*", "", names(seq_ITS2)) %in% gcf_id]
seq_ITS1_filter # 680
seq_ITS2_filter # 688


seq_ITS1_filter_stats <- table(sub("__.*", "", names(seq_ITS1_filter))) %>% as.data.frame() %>% arrange(-Freq)
seq_ITS2_filter_stats <- table(sub("__.*", "", names(seq_ITS2_filter))) %>% as.data.frame() %>% arrange(-Freq)

head(seq_ITS1_filter_stats)
head(seq_ITS2_filter_stats)
nrow(seq_ITS1_filter_stats) # 597
nrow(seq_ITS2_filter_stats) # 605



path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/ITS_CGF"
writeXStringSet(seq_ITS1_filter, file.path(path, "CGF_ITS1.fasta"), format = "fasta")
writeXStringSet(seq_ITS2_filter, file.path(path, "CGF_ITS2.fasta"), format = "fasta")


# run blast_CGF.sh



CGF_ITS1_blast <- read.delim("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/blast/CGF_ITS1_blast.tsv")
CGF_ITS2_blast <- read.delim("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/blast/CGF_ITS2_blast.tsv")


CGF_ITS1_blast_filter <- CGF_ITS1_blast %>%
    mutate(genome = sub("__.*", "", qseqid),
           marker = "ITS1") %>%
    relocate(genome, marker) %>%
    group_by(genome) %>%
    slice_max(order_by = pident, n = 1, with_ties = FALSE) %>%
    ungroup()

CGF_ITS2_blast_filter <- CGF_ITS2_blast %>%
    mutate(genome = sub("__.*", "", qseqid),
           marker = "ITS2") %>%
    relocate(genome, marker) %>%
    group_by(genome) %>%
    slice_max(order_by = pident, n = 1, with_ties = FALSE) %>%
    ungroup()

head(CGF_ITS1_blast_filter)
head(CGF_ITS2_blast_filter)

nrow(CGF_ITS1_blast_filter) # 592
nrow(CGF_ITS2_blast_filter) # 603

table(CGF_ITS1_blast_filter$pident)
table(CGF_ITS2_blast_filter$pident)


GCF_blast <- rbind(CGF_ITS1_blast_filter, CGF_ITS2_blast_filter)
table(GCF_blast$marker)



GCF_blast_plot <- GCF_blast %>%
    select(genome, marker, pident) %>%
    mutate(pident_plot = ifelse(pident < 90, 90, pident)) %>%
    mutate(marker = factor(marker, levels = c("ITS2", "ITS1")))



prop_gcf <- GCF_blast_plot %>%
  group_by(marker) %>%
  summarise(
    proportion = sum(pident >= 98.5) / n() * 100,
    label = paste0("marker: ", unique(marker), " (n = ", n(), ")", "\n", round(proportion, 1), "% seq with pident >= 98.5")
  )
prop_gcf


his_gcf <- ggplot(GCF_blast_plot, aes(x = pident_plot, y = marker, fill = marker)) +
  geom_density_ridges(
    alpha = 0.5,
    scale = 1,
    color = "black",
    rel_min_height = 0,
    bandwidth = 0.15) +
  scale_fill_manual(values = c("ITS1" = "grey60", "ITS2" = "black")) +
  geom_text(
    data = prop_gcf,
    aes(x = 90, y = marker, label = label),
    hjust = 0,
    nudge_y = 0.7,
    size = 3) +
  labs(
    x = "% Identity to the Atlas",
    y = "Density",
    fill = NULL,
    color = NULL) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.key.size = unit(0.15, "in"),
            legend.position = "none")


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/analysis"
ggsave(file.path(path, "his_CGF_blast.pdf"), his_gcf, width = 3, height = 3)








head(GCF_blast)


GCF_blast_filter <- GCF_blast %>% filter(pident >= 98.5)
df_atlas <- data.frame(id = unique(GCF_blast_filter$sseqid))
nrow(df_atlas)


atlas_blast <- read.delim("/data/wangxinyu/Fungi_Atlas/Analysis/Main/blast/fungi/refseq_filterdp_derep_unite_final.tsv")
head(atlas_blast)

df_atlas_meta <- merge(df_atlas, atlas_blast, by.x = "id", by.y = "qseqid", all.x = T)
df_atlas_meta_plot <- df_atlas_meta %>%
    filter(!is.na(pident)) %>%
    mutate(marker = sub(".*(ITS[12]).*", "\\1", id)) %>%
    select(id, marker, pident, qcovs) %>%
    mutate(pident_plot = ifelse(pident < 90, 90, pident),
           qcovs_plot = ifelse(qcovs < 90, 90, qcovs)) %>%
    mutate(marker = factor(marker, levels = c("ITS2", "ITS1")))

head(df_atlas_meta)
nrow(df_atlas_meta) # 334


prop_atlas <- df_atlas_meta_plot %>%
  group_by(marker) %>%
  summarise(
    proportion = sum(pident < 98.5) / n() * 100,
    label = paste0("marker: ", unique(marker)," (n = ", n(), ")", "\n", round(proportion, 1), "% seq with pident < 98.5")
  )

prop_atlas


his_atlas <- ggplot(df_atlas_meta_plot, aes(x = pident_plot, y = marker, fill = marker)) +
  geom_density_ridges(
    alpha = 0.5,
    scale = 1,
    color = "black",
    rel_min_height = 0,
    bandwidth = 0.15) +
  scale_fill_manual(values = c("ITS1" = "grey60", "ITS2" = "black")) +
  geom_text(
    data = prop_atlas,
    aes(x = 90, y = marker, label = label),
    hjust = 0,
    nudge_y = 0.7,
    size = 3) +
  labs(
    x = "% Identity to the UNITE",
    y = "Density",
    fill = NULL,
    color = NULL) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.key.size = unit(0.15, "in"),
            legend.position = "none")


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/analysis"
ggsave(file.path(path, "his_CGF_blast_atlas.pdf"), his_atlas, width = 3, height = 3)




path < "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/analysis"
write.csv(GCF_blast, file.path(path, "GCF_blast.csv"))
write.csv(df_atlas_meta, file.path(path, "GCF_blast_atlas.csv"))
