# Purpose: Compare three workflow outputs with mock ground truth for diversity, detection, abundance, and composition.
# Panels: Figure S2B-E.

library(phyloseq)
library(microbiome)
library(tidyverse)
library(qiime2R)
library(vegan)
library(ggpubr)



color_t <- "black"
color_m <- "#2f3e46"
color_q <- "#cad2c5"
color_d <- "#a5a58d"




tax <- read.csv("/home/wangxy/ITS_benchmark/data/blast/refseq_fun_obipcr_its2_v97_blast_filter.csv", row.names = 1)
tax <- tax %>%
  select(qseqid, sseqid)%>%
  separate(sseqid,
           into = c("species_name", "accession", "sh", "type", "taxonomy"),
           sep = "\\|") %>%
  separate(taxonomy,
           into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";") %>%
  select(-c("species_name", "accession", "sh", "type")) %>%
  dplyr::rename("seq" = "qseqid")


head(tax)



process_mock <- function(file, tax) {

  otu <- read.table(file)

  otu <- otu %>%
    mutate(seq = str_extract(V1, "NR_[^_]+_sub\\[[^]]+\\]")) %>%
    select(-V1) %>%
    dplyr::rename("abundance" = "V2") %>%
    group_by(seq) %>%
    summarise(abundance = sum(abundance), .groups = "drop")


  div <- data.frame(
    Observed = 50,
    Shannon = vegan::diversity(otu$abundance, index = "shannon"),
    Simpson = vegan::diversity(otu$abundance, index = "simpson")
  )


  otu_tax <- merge(otu, tax, by = "seq", all.x = TRUE)

  # species level
  otu_tax_species <- otu_tax %>%
    select(species, abundance) %>%
    group_by(species) %>%
    summarise(abundance = sum(abundance), .groups = "drop")

  # genus level
  otu_tax_genus <- otu_tax %>%
    select(genus, abundance) %>%
    group_by(genus) %>%
    summarise(abundance = sum(abundance), .groups = "drop")


  sample_name <- tools::file_path_sans_ext(basename(file))
  sample_name <- sub("_abundance$", "", sample_name)

  div$sample <- sample_name
  otu_tax_genus$sample <- sample_name

  return(list(div = div, otu = otu_tax_genus))

  cat("done:", sample_name, "\n")
}



files <- list.files("/home/wangxy/ITS_benchmark/data/fq", pattern = "abundance.txt$", full.names = TRUE)
files



div_true <- bind_rows(lapply(files, function(f) process_mock(f, tax)$div))
div_true <- div_true %>% column_to_rownames("sample")
colnames(div_true) <- paste0(colnames(div_true), "_true")
head(div_true)



otu_true <- bind_rows(lapply(files, function(f) process_mock(f, tax)$otu))
otu_true

otu_true <- otu_true %>%
  pivot_wider(
    names_from = genus,
    values_from = abundance,
    values_fill = 0,
    values_fn = sum
  ) %>%
  column_to_rownames("sample") %>%
  as.data.frame()

dim(otu_true)
otu_true[1:4, 1:4]
rowSums(otu_true)






ps_mycogap <- readRDS("/home/wangxy/ITS_benchmark/process/mycogap/3_phyloseq/all/mock_ps_all.rds")
ps_mycogap


div_mycogap <- estimate_richness(ps_mycogap, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
colnames(div_mycogap) <- paste0(colnames(div_mycogap), "_mycogap")
head(div_mycogap)


ps_mycogap <- transform_sample_counts(ps_mycogap, function(x) x / sum(x))
ps_mycogap <- aggregate_taxa(ps_mycogap, level = "Genus")
ps_mycogap

otu_mycogap <- t(as.data.frame(otu_table(ps_mycogap))) %>% as.data.frame()
otu_mycogap[1:4, 1:4]




ps_qiime2 <- qza_to_phyloseq(features = "/home/wangxy/ITS_benchmark/process/qiime2/process/mock_seqtab_cluster.qza",
                             taxonomy = "/home/wangxy/ITS_benchmark/process/qiime2/process/mock_taxtab_cluster.qza")
ps_qiime2


div_qiime2 <- estimate_richness(ps_qiime2, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
colnames(div_qiime2) <- paste0(colnames(div_qiime2), "_qiime2")
head(div_qiime2)


ps_qiime2 <- transform_sample_counts(ps_qiime2, function(x) x / sum(x))
ps_qiime2 <- aggregate_taxa(ps_qiime2, level = "Genus")
ps_qiime2

otu_qiime2 <- t(as.data.frame(otu_table(ps_qiime2))) %>% as.data.frame()
colnames(otu_qiime2) <- paste0("g__", colnames(otu_qiime2))
otu_qiime2[1:4, 1:4]




ps_dada2 <- readRDS("/home/wangxy/ITS_benchmark/process/dada2/mock_dada2_its.rds")
ps_dada2


div_dada2 <- estimate_richness(ps_dada2, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
colnames(div_dada2) <- paste0(colnames(div_dada2), "_dada2")
head(div_dada2)


ps_dada2 <- transform_sample_counts(ps_dada2, function(x) x / sum(x))
ps_dada2 <- aggregate_taxa(ps_dada2, level = "Genus")
ps_dada2

otu_dada2 <- t(as.data.frame(otu_table(ps_dada2))) %>% as.data.frame()
otu_dada2[1:4, 1:4]






div_all <- div_true %>%
  rownames_to_column("sample") %>%
  left_join(rownames_to_column(div_mycogap, "sample"), by = "sample") %>%
  left_join(rownames_to_column(div_qiime2, "sample"), by = "sample") %>%
  left_join(rownames_to_column(div_dada2, "sample"), by = "sample")
head(div_all)

comparisons <- list(
  c("true", "mycogap"),
  c("true", "qiime2"),
  c("true", "dada2"),
  c("mycogap", "qiime2"),
  c("mycogap", "dada2"),
  c("qiime2", "dada2")
)

div_stats <- do.call(rbind, lapply(c("Observed", "Shannon", "Simpson"), function(m) {
  do.call(rbind, lapply(comparisons, function(comp) {
    g1 <- paste0(m, "_", comp[1])
    g2 <- paste0(m, "_", comp[2])

    test <- t.test(div_all[[g1]], div_all[[g2]], paired = TRUE)

    mean1 <- mean(div_all[[g1]], na.rm = TRUE)
    mean2 <- mean(div_all[[g2]], na.rm = TRUE)
    percent_change <- round(((mean2 - mean1) / mean1) * 100, 1)

    data.frame(
      metric = m,
      group1 = comp[1],
      group2 = comp[2],
      mean1 = mean1,
      mean2 = mean2,
      percent_change = percent_change,
      p_value = test$p.value
    )
  }))
}))

div_stats



div_all_change <- div_all %>%
  mutate(
    Observed_per_change_mycogap = (Observed_mycogap - Observed_true) / Observed_true * 100,
    Observed_per_change_qiime2  = (Observed_qiime2 - Observed_true) / Observed_true * 100,
    Observed_per_change_dada2   = (Observed_dada2 - Observed_true) / Observed_true * 100,

    Shannon_per_change_mycogap = (Shannon_mycogap - Shannon_true) / Shannon_true * 100,
    Shannon_per_change_qiime2  = (Shannon_qiime2 - Shannon_true) / Shannon_true * 100,
    Shannon_per_change_dada2   = (Shannon_dada2 - Shannon_true) / Shannon_true * 100,

    Simpson_per_change_mycogap = (Simpson_mycogap - Simpson_true) / Simpson_true * 100,
    Simpson_per_change_qiime2  = (Simpson_qiime2 - Simpson_true) / Simpson_true * 100,
    Simpson_per_change_dada2   = (Simpson_dada2 - Simpson_true) / Simpson_true * 100
  )

head(div_all_change)



path <- "/home/wangxy/ITS_benchmark/analysis"
write.csv(div_all_change, file.path(path, "div.csv"))
write.csv(div_stats, file.path(path, "div_stats.csv"))



div_all_change_long <- div_all_change %>%
 select(sample, ends_with("change_mycogap"), ends_with("change_qiime2"), ends_with("change_dada2")) %>%
  pivot_longer(
    cols = -sample,
    names_to = c("metric", "method"),
    names_pattern = "(.*)_per_change_(.*)",
    values_to = "value") %>%
  mutate(method = factor(method, levels = c("mycogap", "qiime2", "dada2")))


head(div_all_change_long)



box <- ggplot(div_all_change_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.7, color = "black", outlier.shape = 16) +
  scale_x_discrete(labels = c("mycogap" = "MycoGAP", "qiime2" = "QIIME2_ITS", "dada2" = "DADA2_ITS")) +
  scale_y_continuous(labels = function(x) signif(x, 2),
                     expand = expansion(mult = c(0.05, 0.11))) +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw() +
  labs(y = "Relative bias (%)") +
  scale_fill_manual(values = c("mycogap" = color_m, "qiime2" = color_q, "dada2" = color_d)) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = list(c("mycogap", "qiime2"), c("mycogap", "dada2")),
    label = "p.format",
    #vjust = -0.5,
    size = 3,
    tip.length = 0,
    bracket.size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

path <- "/home/wangxy/ITS_benchmark/analysis"
ggsave(file.path(path, "box_div.pdf"), box, width = 3.5, height = 3)




otu_true[1:4, 1:4]
otu_mycogap[1:4, 1:4]
otu_qiime2[1:4, 1:4]
otu_dada2[1:4, 1:4]


otu_all <- bind_rows(
  true = otu_true %>% rownames_to_column("sample"),
  mycogap = otu_mycogap %>% rownames_to_column("sample"),
  qiime2 = otu_qiime2 %>% rownames_to_column("sample"),
  dada2 = otu_dada2 %>% rownames_to_column("sample"),
  .id = "group") %>% arrange(sample, group)

otu_all[is.na(otu_all)] <- 0
otu_all[1:5, 1:5]

path <- "/home/wangxy/ITS_benchmark/analysis"
write.csv(otu_all, file.path(path, "feature_table.csv"))



run_permanova <- function(data, method) {


  otu_sub <- data %>%
    filter(group %in% c("true", method)) %>%
    mutate(id = paste(group, sample, sep = "_")) %>%
    column_to_rownames("id") %>%
    select(-group, -sample)

  # metadata
  meta_sub <- data %>%
    filter(group %in% c("true", method)) %>%
    mutate(
      id = paste(group, sample, sep = "_"),
      group = factor(group, levels = c("true", method))
    ) %>%
    select(id, sample, group)

  # Bray-Curtis distance
  dist_sub <- vegdist(otu_sub, method = "bray")


  meta_sub <- meta_sub %>% slice(match(rownames(as.matrix(dist_sub)), id))
  stopifnot(all(rownames(as.matrix(dist_sub)) == meta_sub$id))

  # PERMANOVA
  set.seed(6311)
  perm_res <- adonis2(dist_sub ~ group, data = meta_sub, permutations = 999)

  return(perm_res)
}

# run
perm_m <- run_permanova(otu_all, "mycogap")
perm_q <- run_permanova(otu_all, "qiime2")
perm_d <- run_permanova(otu_all, "dada2")
perm_m
perm_q
perm_d


otu_all2 <- otu_all %>%
    mutate(id = paste(group, sample, sep = "_")) %>%
    column_to_rownames("id") %>%
    select(-group, -sample)

otu_all2[1:5, 1:5]


dist_bc <- vegdist(otu_all2, method = "bray")
str(dist_bc)

# PCoA
pcoa_res <- cmdscale(dist_bc, k = 2, eig = TRUE)


pcoa_df <- data.frame(
  PCoA1 = pcoa_res$points[, 1],
  PCoA2 = pcoa_res$points[, 2]
) %>%
  rownames_to_column("id") %>%
  separate(id, into = c("group", "sample1", "sample2"), sep = "_") %>%
  mutate(sample = paste(sample1, sample2, sep = "_")) %>%
  select(PCoA1, PCoA2, group, sample)

pcoa_df

eig <- pcoa_res$eig
eig_pos <- eig[eig > 0]
var_exp <- round(eig_pos[1:2] / sum(eig_pos) * 100, 1)
var_exp



anno_label <- paste0(
  "PERMANOVA\n",
  "MycoGAP: p-value = ", signif(perm_m$`Pr(>F)`[1], 3), "\n",
  "QIIME2_ITS: p-value = ", signif(perm_q$`Pr(>F)`[1], 3), "\n",
  "DADA2_ITS: p-value = ", signif(perm_d$`Pr(>F)`[1], 3)
)
anno_label

anno_label <- "All pairwise PERMANOVA vs. true data: p = 1"
anno_label <- ""


pcoa <- ggplot() +
  geom_point(data = pcoa_df %>% filter(group != "true"),
             aes(x = PCoA1, y = PCoA2, color = group), size = 1, shape = 16, alpha = 0.8) +
  geom_point(data = pcoa_df %>% filter(group == "true"),
             aes(x = PCoA1, y = PCoA2), color = "#222222", shape = 17, size = 1, alpha = 0.8,
             show.legend = FALSE) +
  scale_color_manual(values = c("mycogap" = color_m, "qiime2" = color_q, "dada2" = color_d),
                     labels = c("mycogap" = "MycoGAP", "qiime2" = "QIIME2_ITS", "dada2" = "DADA2_ITS"),
                     breaks = c("mycogap", "qiime2", "dada2")) +
  labs(x = paste0("PCoA1 (", var_exp[1], "%)"),
       y = paste0("PCoA2 (", var_exp[2], "%)"),
       color = NULL) +
  guides(color = guide_legend(ncol = 3, override.aes = list(size = 2))) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = "bottom",
        legend.key.size = unit(0.15, "in"),
        legend.key.width = unit(0.05, "in"))


path <- "/home/wangxy/ITS_benchmark/analysis"
ggsave(file.path(path, "pcoa.pdf"), pcoa, width = 3, height = 3)






dist_res <- as.matrix(dist_bc) %>%
  as_tibble(rownames = "s1") %>%
  pivot_longer(-s1, names_to = "s2", values_to = "distance") %>%
  mutate(
    method1 = str_extract(s1, "^[^_]+"),
    mock1   = str_extract(s1, "mock_\\d+"),
    method2 = str_extract(s2, "^[^_]+"),
    mock2   = str_extract(s2, "mock_\\d+")
  ) %>%
  filter(mock1 == mock2,
         method2 == "true",
         method1 %in% c("mycogap", "qiime2", "dada2")) %>%
  select(sample = mock1, method = method1, distance) %>%
  pivot_wider(names_from  = method,
              values_from = distance,
              names_prefix = "dist_") %>%
  arrange(sample)

dist_res



cor_res <- otu_all %>%
  group_by(sample) %>%
  group_modify(~ {

    df_sample <- .x %>%
      column_to_rownames("group")

    truth_abund   <- as.numeric(df_sample["true", ])
    mycogap_abund <- as.numeric(df_sample["mycogap", ])
    qiime2_abund  <- as.numeric(df_sample["qiime2", ])
    dada2_abund   <- as.numeric(df_sample["dada2", ])

    valid_idx_mycogap <- (truth_abund + mycogap_abund) > 0
    valid_idx_qiime2  <- (truth_abund + qiime2_abund)  > 0
    valid_idx_dada2   <- (truth_abund + dada2_abund)   > 0

    # mycogap vs true
    test_m <- cor.test(
      truth_abund[valid_idx_mycogap],
      mycogap_abund[valid_idx_mycogap],
      method = "spearman"
    )

    # qiime2 vs true
    test_q <- cor.test(
      truth_abund[valid_idx_qiime2],
      qiime2_abund[valid_idx_qiime2],
      method = "spearman"
    )

    # dada2 vs true
    test_d <- cor.test(
      truth_abund[valid_idx_dada2],
      dada2_abund[valid_idx_dada2],
      method = "spearman"
    )

    tibble(
      cor_mycogap = unname(test_m$estimate),
      #p_mycogap = test_m$p.value,

      cor_qiime2 = unname(test_q$estimate),
      #p_qiime2 = test_q$p.value

      cor_dada2 = unname(test_d$estimate),
      #p_dada2 = test_d$p.value
    )
  }) %>%
  ungroup()

cor_res




com_res <- merge(cor_res, dist_res, by = "sample")
com_res

path <- "/home/wangxy/ITS_benchmark/analysis"
write.csv(com_res, file.path(path, "com_stats.csv"))


com_long <- com_res %>%
  pivot_longer(
    cols = -sample,
    names_to = "metric_method",
    values_to = "value"
  ) %>%
  separate(metric_method, into = c("metric", "group"), sep = "_") %>%
  mutate(
    metric = factor(
      metric,
      levels = c("cor", "dist"),
      labels = c("Spearman rho", "Bray-Curtis")
    )
  ) %>%
  mutate(group = factor(group, levels = c("mycogap", "qiime2", "dada2")))

head(com_long)



box <- ggplot(com_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(alpha = 0.7, color = "black", outlier.shape = 16) +
  facet_wrap(~metric, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.11))) +
  scale_x_discrete(labels = c("mycogap" = "MycoGAP", "qiime2" = "QIIME2_ITS", "dada2" = "DADA2_ITS")) +
  theme_bw() +
  labs(y = NULL) +
  scale_fill_manual(values = c("true" = color_t, "mycogap" = color_m, "qiime2" = color_q, "dada2" = color_d)) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = list(c("mycogap", "qiime2"), c("mycogap", "dada2")),
    label = "p.format",
    #vjust = -0.5,
    size = 3,
    tip.length = 0,
    bracket.size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

path <- "/home/wangxy/ITS_benchmark/analysis"
ggsave(file.path(path, "box_com.pdf"), box, width = 2.5, height = 3)





otu_all[1:5, 1:5]


otu_all_long <- otu_all %>%
  pivot_longer(
    -c(group, sample),
    names_to = "taxon",
    values_to = "abund"
  ) %>%
  mutate(present = abund > 0)
head(otu_all_long)


truth <- otu_all_long %>%
  filter(group == "true") %>%
  select(sample, taxon, true_present = present, true_abund = abund)
head(truth)


pred <- otu_all_long %>%
  filter(group != "true") %>%
  rename(method = group, pred_present = present, pred_abund = abund)
head(pred)


taxa_res <- pred %>%
  left_join(truth, by = c("sample", "taxon")) %>%
  group_by(method, sample) %>%
  summarise(

    TP = sum(pred_present & true_present),
    FP = sum(pred_present & !true_present),
    FN = sum(!pred_present & true_present),


    Precision = TP / (TP + FP),
    Recall = TP / (TP + FN),
    F1_score = 2 * Precision * Recall / (Precision + Recall),

    .groups = "drop"
  )

taxa_res


path <- "/home/wangxy/ITS_benchmark/analysis"
write.csv(taxa_res, file.path(path, "taxa_stats.csv"))



taxa_res_long <- taxa_res %>%
  select(method, sample, Precision, Recall, F1_score) %>%
  pivot_longer(
    cols = c(Precision, Recall, F1_score),
    names_to = "metric",
    values_to = "value"
  ) %>%

  mutate(
    metric = factor(metric, levels = c("Precision", "Recall", "F1_score")),
    method = factor(method, levels = c("mycogap", "qiime2", "dada2"))
  )

head(taxa_res_long)


box <- ggplot(taxa_res_long, aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.7, color = "black", outlier.shape = 16) +
  facet_wrap(~metric) +
  #scale_y_continuous(limits = c(0.75, 1.02)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.11))) +
  scale_x_discrete(labels = c("mycogap" = "MycoGAP", "qiime2" = "QIIME2_ITS", "dada2" = "DADA2_ITS")) +
  theme_bw() +
  labs(y = NULL) +
  scale_fill_manual(values = c("true" = color_t, "mycogap" = color_m, "qiime2" = color_q, "dada2" = color_d)) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = list(c("mycogap", "qiime2"), c("mycogap", "dada2")),
    label = "p.format",
    #vjust = -0.5,
    size = 3,
    tip.length = 0,
    bracket.size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

path <- "/home/wangxy/ITS_benchmark/analysis"
ggsave(file.path(path, "box_taxa.pdf"), box, width = 3, height = 3)





# - more detailed check -
taxa_stats <- read.csv("/home/wangxy/ITS_benchmark/analysis/taxa_stats.csv", row.names = 1)
head(taxa_stats)

taxa_stats_summary <- taxa_stats %>%
  group_by(method) %>%
  summarise(
    mean_Precision = mean(Precision, na.rm = TRUE),
    mean_Recall = mean(Recall, na.rm = TRUE),
    mean_F1_score = mean(F1_score, na.rm = TRUE)
  )

taxa_stats_summary


com_stats <- read.csv("/home/wangxy/ITS_benchmark/analysis/com_stats.csv", row.names = 1)
head(com_stats)

round(mean(com_stats$cor_mycogap), 3) # 0.895
round(mean(com_stats$cor_qiime2), 3) # 0.785
round(mean(com_stats$cor_dada2), 3) # 0.893

round(mean(com_stats$dist_mycogap), 3) # 0.068
round(mean(com_stats$dist_qiime2), 3) # 0.122
round(mean(com_stats$dist_dada2), 3) # 0.07
