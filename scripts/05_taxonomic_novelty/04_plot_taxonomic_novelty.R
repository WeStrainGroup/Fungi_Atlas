# Purpose: Plot database-match identity distributions and phylum-stratified operational novelty levels.
# Panels: Figure 2E-F.

library(tidyverse)
library(phyloseq)
library(microbiome)
library(data.table)
library(vegan)
library(ggrepel)
library(ggpubr)
library(scales)
library(gridExtra)
library(paletteer)
library(igraph)
library(ggraph)
library(tidygraph)
library(ComplexUpset)
library(eulerr)
library(patchwork)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(ggbeeswarm)
library(broom)
library(ggridges)

# --- Visualize BLAST results for ASVs ---
# Unite
blast_ITS <- read.delim("/data/wangxinyu/ITS_Public/Project/Analysis/blast/fungi/refseq_filterdp_derep_unite_final.tsv", header = TRUE)

# clean
blast_ITS <- blast_ITS %>% select (qseqid, pident)
blast_ITS$category <- "ITS"
head(blast_ITS)
nrow(blast_ITS) # 55902

# proportion of novel taxa
per1_ITS <- round(mean(blast_ITS$pident < 98.5) * 100, 1)
per1_ITS

# load total read counts per ASV
ASV_counts0_ITS <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/refseq_filterdp_derep_counts.csv", row.names = 1)
head(ASV_counts0_ITS)
nrow(ASV_counts0_ITS) # 57197

# merge with BLAST results
blast_ITS$qseqid <- gsub(";size=\\d+$", "", blast_ITS$qseqid) # remove qseqid suffix
head(blast_ITS)
nrow(blast_ITS)

ASV_counts_ITS <- merge(blast_ITS, ASV_counts0_ITS, by.x = "qseqid", by.y = "row.names", all.x = TRUE)
head(ASV_counts_ITS)


# proportion of reads from novel taxa
count2_ITS <- ASV_counts_ITS %>%
  filter(pident < 98.5) %>%
  summarise(sum_counts = sum(counts, na.rm = TRUE)) %>%
  pull(sum_counts)
count2_ITS

per2_ITS <- round(count2_ITS / sum(ASV_counts_ITS$counts) * 100, 1)
per2_ITS

# 16S
# Silva
blast_16S <- read.delim("/data/wangxinyu/ITS_Public/Project/Analysis/blast/cell2025/analysis/cell2025_asv_filter2_0.05_silva_final.tsv", header = TRUE)

# clean
blast_16S <- blast_16S %>% select (qseqid, pident)
blast_16S$category <- "16S"
head(blast_16S)
nrow(blast_16S)

# proportion of novel ASVs
per1_16S <- round(mean(blast_16S$pident < 98.5) * 100, 1)
per1_16S

# compute read proportions by identity threshold
# load total read counts per ASV
ASV_counts0_16S <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/blast/cell2025/analysis/cell2025_asv_filter_counts.csv", row.names = 1)
head(ASV_counts0_16S)
nrow(ASV_counts0_16S)

# merge with BLAST results
ASV_counts_16S <- merge(blast_16S, ASV_counts0_16S, by.x = "qseqid", by.y = "ASV", all.x = TRUE)
head(ASV_counts_16S)
nrow(ASV_counts_16S)


# proportion of reads from novel taxa
count2_16S <- ASV_counts_16S %>%
  filter(pident < 98.5) %>%
  summarise(sum_counts = sum(counts, na.rm = TRUE)) %>%
  pull(sum_counts)

count2_16S

per2_16S <- round(count2_16S / sum(ASV_counts_16S$counts) * 100, 1)
per2_16S


# merge datasets
blast_merge <- rbind(blast_ITS, blast_16S)
head(blast_merge)
tail(blast_merge)

# colors
col_ITS <- "#1e1e24"
col_16S <- "#8ab0ab"

col <- c("ITS" = col_ITS,
         "16S" = col_16S)

# visualization
# build annotation text table
table <- paste0("                             ", "ITS", "  ", "16S",
                "\n% of novel ASVs: ", per1_ITS, "  ", per1_16S,
                "\n% of novel reads: ", per2_ITS, "  ", per2_16S)
table

# visualization
his_blast <- ggplot(blast_merge, aes(x = pident, fill = category)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, position = "identity", alpha = 0.5, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = col, breaks = c("ITS", "16S")) + # fixed order
  coord_cartesian(xlim = c(72.5, NA)) +
  scale_x_continuous(breaks = seq(70, 100, by = 5)) +
      labs(x = "Percent identity",
           y = "Density",
           color = NULL) +
      geom_vline(xintercept = 98.5, size = 0.5, linetype = "dashed") +
      annotate("text", x = 75, y = 0.7,
           label = table,
           size = 3, hjust = 0) +
      annotate("text", x = 97.5, y = 0.4,
           label = paste0("threshold: 98.5"),
           size = 3, hjust = 1) +
      theme_bw() +
      theme(panel.border = element_rect(color = "black", linewidth = 1),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.title = element_blank(),
            legend.key.size = unit(0.15, "in"),
            legend.position = c(0.2, 0.2))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "his_blast_unite_silva.pdf"), his_blast, width = 3, height = 3)



# --- Novelty of newly observed ASVs ---
df_level <- read.csv("/data/huangkailang/project/fungi_ITS/data_for_plot/new_ASV/ASV_number_annotation_level.csv") %>% arrange(-Freq)
head(df_level)
unique(df_level$phylum)

df_level <- df_level %>% filter(phylum != "Fungi_phy_Incertae_sedis")


df_level <- df_level %>%
  mutate(Novelty_level = case_when(
    Annotation_level == "Species" ~ "Known species",
    Annotation_level == "Genus" ~ "Novel species",
    Annotation_level == "Family" ~ "Novel genus",
    Annotation_level == "Order" ~ "Novel family",
    Annotation_level == "Class" ~ "Novel order",
    Annotation_level == "Phylum" ~ "Novel class",
  ))

df_level$Novelty_level <- factor(df_level$Novelty_level,
                                 levels = c("Known species", "Novel species", "Novel genus", "Novel family", "Novel order", "Novel class"))
df_level$phylum <- factor(df_level$phylum, levels = unique(df_level$phylum))

p <-
  ggplot(df_level, aes(x = phylum, y = Novelty_level, size = Freq)) + #fill = log10(Freq + 1),
  geom_point(alpha = 0.7, shape = 16, fill = "#1e1e24") +
  labs(x = "Assigned phylum") +
  scale_size(name = "# ASVs",
             range = c(1, 7),
             breaks = c(200, 2000, 20000),
             labels = c(200, 2000, 20000)) +
  theme_bw()+
  theme(legend.position = c(0.735, 0.835),
        legend.direction = "horizontal",
        legend.key.size = unit(0.15, "in"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  guides(size = guide_legend(title.position = "left", label.position = "right"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bubble_ASV_assign.pdf"), p, width = 6, height = 3)
