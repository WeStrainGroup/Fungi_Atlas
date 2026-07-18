# Purpose: Plot phylum/genus composition, cross-site MDS, and read-based rarefaction/extrapolation.
# Panels: Figure 2A-B; Figure S5A-C.

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

# ---- Composition across all samples at the Phylum level ----
ps <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/ps_merged_ASV_filterdp.rds")
ps_phy_0.01 <- aggregate_rare(ps, level = "Phylum", detection = 0, prevalence = 0.01)
unique(rownames(tax_table(ps_phy_0.01)))

ps_phy_0.01_com <- microbiome::transform(ps_phy_0.01, "compositional")
otu_phy <- as.data.frame(otu_table(ps_phy_0.01_com))
rownames(otu_phy) <- gsub("^p__", "", rownames(otu_phy))
otu_phy[1:5, 1:5]
ncol(otu_phy) # 37417

# rename 'Other'
rownames(otu_phy)[rownames(otu_phy) == "Other"] <- "Other (prev. < 1%)"
write.csv(otu_phy, "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/ra_phy_0.01.csv")
otu_phy <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/ra_phy_0.01.csv", row.names = 1)


# reorder samples by Ascomycota abundance
ascomycota_abundance <- as.numeric(otu_phy["Ascomycota", ])
otu_phy <- otu_phy[, order(-ascomycota_abundance)]

# convert to long format for ggplot
otu_long <- reshape2::melt(as.matrix(otu_phy),
                 varnames = c("phylum", "sample"),
                 value.name = "abundance")

# order Phyla by total abundance (desc)
phylum_total <- tapply(otu_long$abundance, otu_long$phylum, sum)
phylum_order <- names(sort(phylum_total, decreasing = TRUE))
otu_long$phylum <- factor(otu_long$phylum, levels = phylum_order)
head(otu_long)

# plot
bar_phy <- ggplot(otu_long, aes(x = sample, y = abundance, fill = phylum)) +
  geom_bar(stat = "identity", width = 1, alpha = 0.85) +
  scale_fill_paletteer_d("ggthemes::stata_s2color") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 2)) +  # legend in 2 columns
  theme_bw() +
  labs(x = "Samples (n = 37,417, clean reads of microfungi >= 10,000)",
       y = "Relative abundance",
       fill = "phylum (ordered by abundance)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.3, 0.76),
        legend.key.size = unit(0.15, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_phy.png"), bar_phy, width = 6, height = 3, dpi = 900)
#ggsave(file.path(path, "bar_phy.pdf"), bar_phy, width = 6, height = 3)





# ---- Composition across all samples at the Genus level ----
ps <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/ps_merged_ASV_filterdp.rds")
ps_gen_0.25 <- aggregate_rare(ps, level = "Genus", detection = 0, prevalence = 0.25)
unique(rownames(tax_table(ps_gen_0.25))) # 19

ps_gen_0.25_com <- microbiome::transform(ps_gen_0.25, "compositional")
otu_gen <- as.data.frame(otu_table(ps_gen_0.25_com))
rownames(otu_gen) <- gsub("^g__", "", rownames(otu_gen))
otu_gen[1:5, 1:5]
ncol(otu_gen) # 37417

# rename 'Other' and 'Unknown'
rownames(otu_gen)[rownames(otu_gen) == "Other"] <- "Other (prev. < 25%)"
rownames(otu_gen)[rownames(otu_gen) == "Unknown"] <- "Unassigend"

write.csv(otu_gen, "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/ra_gen_0.25.csv")
otu_gen <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/ra_gen_0.25.csv", row.names = 1)
otu_gen[1:3, 1:3]

# reorder samples by Saccharomyces abundance
Saccharomyces_abundance <- as.numeric(otu_gen["Saccharomyces", ])
otu_gen <- otu_gen[, order(-Saccharomyces_abundance)]

# convert to long format for ggplot
otu_long <- reshape2::melt(as.matrix(otu_gen),
                 varnames = c("genus", "sample"),
                 value.name = "abundance")

# order Genera by total abundance (desc)
genus_total <- tapply(otu_long$abundance, otu_long$genus, sum)
genus_order <- names(sort(genus_total, decreasing = TRUE))
otu_long$genus <- factor(otu_long$genus, levels = genus_order)
head(otu_long)

# plot
bar_gen <- ggplot(otu_long, aes(x = sample, y = abundance, fill = genus)) +
  geom_bar(stat = "identity", width = 1, alpha = 0.85) +
  scale_fill_paletteer_d("ggsci::default_igv") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 2)) +  # legend in 2 columns
  theme_bw() +
  labs(x = "Samples (n = 37,417, clean reads of microfungi >= 10,000)",
       y = "Relative abundance",
       fill = "Genus (ordered by abundance)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.15, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_gen.png"), bar_gen, width = 12, height = 3, dpi = 900)












# --- Fungal distribution across sites ---
MDS_point <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged/MDS_site_point.csv", row.names = 1)
head(MDS_point)

MDS_eigen <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged/MDS_site_eigen.csv", row.names = 1)
head(MDS_eigen)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged/meta_gen_p0_filterdp_merged.csv",  row.names = 1)
head(meta)

MDS_point_meta <- MDS_point %>%
  rownames_to_column("id") %>%
  left_join(meta, by = "id")

head(MDS_point_meta)
table(MDS_point_meta$site)

# explained variance
eigen1 <- round(100 * MDS_eigen$eigen[1], 1)
eigen2 <- round(100 * MDS_eigen$eigen[2], 1)
eigen3 <- round(100 * MDS_eigen$eigen[3], 1)
eigen4 <- round(100 * MDS_eigen$eigen[4], 1)


# colors
col <- c("human_gut" = "#457b9d",
         "mouse_gut"  = "#a2142f",
         "environment" = "#edb120")


# MSD1 and 2
MDS <- ggplot() +
  # plot human and environment first
  geom_point(
    data = subset(MDS_point_meta, site != "mouse_gut"),
    aes(x = MDS1, y = MDS2, color = site),
    size = 0.01, shape = 16, alpha = 0.75
  ) +
  # plot mouse_gut last
  geom_point(
    data = subset(MDS_point_meta, site == "mouse_gut"),
    aes(x = MDS1, y = MDS2, color = site),
    size = 0.01, shape = 16, alpha = 0.75
  ) +
  scale_color_manual(values = col,
                    breaks = c("human_gut", "mouse_gut", "environment")) +
  geom_point(size = 0.01, shape = 16, alpha = 0.75) +
  labs(x = paste0("MDS1 (", eigen1, "%)"),
       y = paste0("MDS2 (", eigen2, "%)"),
       color = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  annotate("text", x = 3000, y = 350,
           label = paste0("ANOSIM statistic R: 0.46", "\nsignificance: 0.001"),
           size = 3, hjust = 1) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid = element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.75, 0.15),
        legend.key.size = unit(0.15, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "MDS_site1.pdf"), MDS, width = 3, height = 3)


# MSD3 and 4
MDS <- ggplot() +
  # plot human and environment first
  geom_point(
    data = subset(MDS_point_meta, site != "mouse_gut"),
    aes(x = MDS3, y = MDS4, color = site),
    size = 0.01, shape = 16, alpha = 0.75
  ) +
  # plot mouse_gut last
  geom_point(
    data = subset(MDS_point_meta, site == "mouse_gut"),
    aes(x = MDS3, y = MDS4, color = site),
    size = 0.01, shape = 16, alpha = 0.75
  ) +
  scale_color_manual(values = col,
                    breaks = c("human_gut", "mouse_gut", "environment")) +
  geom_point(size = 0.01, shape = 16, alpha = 0.75) +
  labs(x = paste0("MDS3 (", eigen3, "%)"),
       y = paste0("MDS4 (", eigen4, "%)"),
       color = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  annotate("text", x = 500, y = 200,
           label = paste0("ANOSIM statistic R: 0.46", "\nsignificance: 0.001"),
           size = 3, hjust = 1) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        legend.position = c(0.75, 0.15),
        legend.key.size = unit(0.15, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "MDS_site2.pdf"), MDS, width = 3, height = 3)





# --- Reads-based rarefaction visualization ---
# type 1 size-based
df <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/rarefaction/rarefaction_reads_type1.csv", row.names = 1)
head(df)
table(df$Order.q)

# Draw R/E curves by yourself
df.line <- df %>% filter(Method != "Observed")
df.line <- df.line %>%
  mutate(Method = factor(Method, levels = c("Rarefaction", "Extrapolation")))
df$Order.q.label <- factor(
  paste0("q = ", df$Order.q),
  levels = paste0("q = ", sort(unique(df$Order.q)))
)
df.line$Order.q.label <- factor(
  paste0("q = ", df.line$Order.q),
  levels = paste0("q = ", sort(unique(df$Order.q)))
)

p <- ggplot(df, aes(x = x, y = y, group = interaction(Assemblage, Method, Order.q.label))) +
  geom_line(aes(linetype = Method), data = df.line, lwd = 0.1, alpha = 0.5, color = "grey") +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), data = df, fill = "grey", alpha = 0.2) +
  geom_vline(xintercept = 10000, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~ Order.q.label, nrow = 1, scales = "free_y") +
  labs(x = "Number of reads", y = "Hill numbers based on ASVs") +
  theme_bw() +
  theme(legend.position = "none")

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "line_raref_reads_type1.pdf"), p, width = 9, height = 3)



# type2 sample coverage
df <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/rarefaction/rarefaction_reads_type2.csv", row.names = 1)
head(df)

# Draw R/E curves by yourself
df.line <- df %>% filter(Method != "Observed")
df.line <- df.line %>%
  mutate(Method = factor(Method, levels = c("Rarefaction", "Extrapolation")))


p <- ggplot(df, aes(x = x, y = y, group = interaction(Assemblage, Method))) +
  geom_line(aes(linetype = Method), data = df.line, linewidth = 0.1, alpha = 0.5, color = "grey") +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr), data = df, fill = "grey", alpha = 0.2) +
  geom_vline(xintercept = 10000, linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(x = "Number of reads", y = "Sample coverage") +
  theme_bw() +
  theme(legend.position = "none")

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "line_raref_reads_type2.pdf"), p, width = 3, height = 3)
