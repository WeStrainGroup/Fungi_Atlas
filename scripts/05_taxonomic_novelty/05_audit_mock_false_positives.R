# Purpose: Export recovered mock features, combine matched BLAST results, and summarize false-positive counts and identities.
# Panels: Figure S5E-F.


library(phyloseq)
library(Biostrings)
library(tidyverse)
library(ggpubr)



ps <- readRDS("/home/wangxy/ITS_benchmark/process/mycogap/3_phyloseq/all/mock_ps_all.rds")
ps


# get otu
otu <- as(otu_table(ps), "matrix")
otu <- t(otu)
otu[1:3, 1:3]

# get seq
seqs <- refseq(ps)
seqs

out_dir <- "/home/wangxy/ITS_fp_check/data/feature_mycogap"
dir.create(out_dir, showWarnings = FALSE)


for (s in sample_names(ps)) {
  present <- rownames(otu)[otu[, s] > 0]

  sub_seqs <- seqs[present]


  names(sub_seqs) <- paste0(present, ";count=", otu[present, s])

  writeXStringSet(sub_seqs,
                  filepath = file.path(out_dir, paste0(s, "_sequence.fa")))
}


# run blast_to_mock



files <- list.files(
  "/home/wangxy/ITS_fp_check/data/blast_res",
  pattern = "\\.tsv$",
  full.names = TRUE
)

files

blast_res <- map_dfr(files, function(f) {

  prefix <- basename(f) %>% str_remove("_sequence_blast\\.tsv$")
  read_tsv(f, show_col_types = FALSE) %>% mutate(sample = prefix, .before = 1) %>% as.data.frame()
})


head(blast_res)
table(blast_res$pident)




df <- as.data.frame(table(blast_res$sample)) %>% rename("sample" = "Var1", "mycogap_seqs_hit" = "Freq")


meta1 <- read.table("/home/wangxy/ITS_fp_check/data/feature_mycogap/fa_mycogap_stats.tsv", header = T)
meta1 <- meta1 %>%
    mutate(sample = str_remove(file, "_sequence\\.fa$")) %>%
    select(sample, num_seqs) %>%
    rename("mycogap_seqs" = "num_seqs")

head(meta1)


meta2 <- read.table("/home/wangxy/ITS_fp_check/data/feature_derep/fa_derep_stats.tsv", header = T)
meta2 <- meta2 %>%
    mutate(sample = str_remove(file, "_sequence_derep\\.fa$")) %>%
    select(sample, num_seqs) %>%
    rename("mock_seqs_mutation" = "num_seqs")

head(meta2)


meta <- merge(meta1, meta2, by = "sample") %>%
    mutate(mock_seqs = 50)  %>%
    select(sample, mock_seqs, mock_seqs_mutation, mycogap_seqs)
head(meta)

blast_res_stats <- merge(meta, df, by = "sample")
head(blast_res_stats)

sum(blast_res_stats$mycogap_seqs) # 2775
sum(blast_res_stats$mycogap_seqs_hit) # 2772
ratio_hit <- round(sum(blast_res_stats$mycogap_seqs_hit) / sum(blast_res_stats$mycogap_seqs) * 100, 1)
ratio_hit
# 99.9





blast_res_stats_long <- blast_res_stats %>%
  pivot_longer(
    cols = -sample,
    names_to = "metric",
    values_to = "value")

head(blast_res_stats_long)



box <- ggplot(blast_res_stats_long, aes(x = metric, y = value, fill = metric)) +
  geom_boxplot(alpha = 0.6, color = "black", outlier.shape = 16) +
  scale_x_discrete(labels = c("mock_seqs" = "Mock",
                              "mock_seqs_mutation" = "Mock_mutation",
                              "mycogap_seqs" = "MycoGAP",
                              "mycogap_seqs_hit" = "MycoGAP_hit")) +
  scale_y_continuous(labels = function(x) signif(x, 2),
                     expand = expansion(mult = c(0.05, 0.11))) +
  theme_bw() +
  labs(y = "Number of features") +
  scale_fill_manual(values = c("mock_seqs" = "black",
                               "mock_seqs_mutation" = "black",
                               "mycogap_seqs" = "#2f3e46",
                               "mycogap_seqs_hit" = "#2f3e46")) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = list(c("mock_seqs", "mock_seqs_mutation"),
                       c("mock_seqs", "mycogap_seqs"),
                       c("mock_seqs_mutation", "mycogap_seqs"),
                       c("mycogap_seqs", "mycogap_seqs_hit")),
    label = "p.format",
    #vjust = -0.5,
    size = 3,
    tip.length = 0,
    bracket.size = 0.5) +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")




path <- "/home/wangxy/ITS_fp_check/analysis"
ggsave(file.path(path, "box_mock_fp_check.pdf"), box, width = 3, height = 3)




prop <- blast_res %>%
  summarise(
    prop_100 = sum(pident >= 100) / n() * 100,
    prop_99 = sum(pident >= 99) / n() * 100,
    prop_98.5 = sum(pident >= 98.5) / n() * 100,
    label = paste0("n = ", n(), " (", ratio_hit, "% of total features)", "\n",
                round(prop_98.5, 1), "% with pident >= 98.5", "\n",
                round(prop_99, 1), "% with pident >= 99", "\n",
                round(prop_100, 1), "% with pident = 100")
  )
prop



his <- ggplot(blast_res, aes(x = pident)) +
  geom_histogram(aes(y = after_stat(density)), bins = 25, position = "identity", alpha = 0.7, color = "black", linewidth = 0.4) +
  labs(x = "% Identity to features in mock",
       y = "Density",
       color = NULL) +
  annotate("text", x = 98.5, y = 11, label = prop$label,  hjust = 0, size = 3) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())



path <- "/home/wangxy/ITS_fp_check/analysis"
ggsave(file.path(path, "his_mock_fp_check.pdf"), his, width = 3, height = 3)










path <- "/home/wangxy/ITS_fp_check/analysis"
write.csv(blast_res, file.path(path, "blast_res_mock_fp_check.csv"))
write.csv(blast_res_stats, file.path(path, "blast_res_stats_mock_fp_check.csv"))

head(blast_res_stats)
round(mean(blast_res_stats$mock_seqs_mutation)) # 85
round(mean(blast_res_stats$mycogap_seqs)) # 56
