# Purpose: Plot atlas read composition, depth, diversity, prevalence, abundance, and sample-based rarefaction summaries.
# Panels: Figure 1D-F, Figure 2C-D, Figure 2G, and Figure S5D.

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

# --- Visualize read depth distribution ---
depth <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_merged.csv", row.names = 1)
depth <- depth %>% filter(depth_all >= 10000)
head(depth)

# compute medians
median_raw <- round(median(depth$depth_raw), 0)
median_all <- round(median(depth$depth_all), 0)
median_microfungi <- round(median(depth$depth_microfungi),0)
median_raw
median_all
median_microfungi


# visualize all species vs. microfungi
# convert to long format
depth_long <- pivot_longer(
  depth,
  cols = c(depth_all, depth_microfungi),
  names_to = "category",
  values_to = "depth")
head(depth_long)


# replace zeros with ones
depth_long <- depth_long %>%
  mutate(depth = ifelse(depth == 0, 1, depth))

col_all <- "#fcd5ce"
col_microfungi <- "#1e1e24"

col <- c("depth_all" = col_all,
         "depth_microfungi" = col_microfungi)

label <- c("depth_all" = "all species",
           "depth_microfungi" = "microfungi")

his_depth1 <- ggplot(depth_long, aes(x = depth, fill = category)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.5, color = "black", linewidth = 0.25) +
      scale_fill_manual(values = col,
                        labels = label,
                        breaks = c("depth_all", "depth_microfungi")) + # fixed order
      scale_x_log10() +
      coord_cartesian(xlim = c(1, 1e7)) +
      labs(x = "Number of reads",
           y = "Number of samples",
           fill = NULL) +
      theme_bw() +
      # median line
      geom_segment(aes(x = median_all, xend = median_all, y = 0, yend = 15000),
             color = "#ad6a6c", linewidth = 0.5) +
      geom_segment(aes(x = median_microfungi, xend = median_microfungi, y = 0, yend = 10000),
             color = col_microfungi, linewidth = 0.5) +
      # median annotation
      annotate("text", x = 2e1, y = 14000,
           label = paste0("median: ", format(median_all, big.mark = ",")),
           color = "#ad6a6c", size = 3, hjust = 0) +
      annotate("text", x = 1.5e1, y = 9000,
           label = paste0("median: ", format(median_microfungi, big.mark = ",")),
           color = col_microfungi, size = 3, hjust = 0) +
      theme(panel.border = element_rect(color = "black", linewidth = 1),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.position = c(0.25, 0.25),
            legend.key.size = unit(0.15, "in"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "his_depth1.pdf"), his_depth1, width = 3, height = 3)

# visualize raw fastq vs. microfungi
# convert to long format
depth_long <- pivot_longer(
  depth,
  cols = c(depth_raw, depth_microfungi),
  names_to = "category",
  values_to = "depth")
head(depth_long)


# replace zeros with ones
depth_long <- depth_long %>%
  mutate(depth = ifelse(depth == 0, 1, depth))

col_raw <- "#003049"
col_microfungi <- "#1e1e24"

col <- c("depth_raw" = col_raw,
         "depth_microfungi" = col_microfungi)

label <- c("depth_raw" = "raw fastq",
           "depth_microfungi" = "microfungi")

his_depth2 <- ggplot(depth_long, aes(x = depth, fill = category)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.5, color = "black", linewidth = 0.25) +
      scale_fill_manual(values = col,
                        labels = label,
                        breaks = c("depth_raw", "depth_microfungi")) + # fixed order
      scale_x_log10() +
      coord_cartesian(xlim = c(1, 1e7)) +
      labs(x = "Number of reads",
           y = "Number of samples",
           fill = NULL) +
      theme_bw() +
      # median line
      geom_segment(aes(x = median_raw, xend = median_raw, y = 0, yend = 17000),
             color = col_raw, linewidth = 0.5) +
      geom_segment(aes(x = median_microfungi, xend = median_microfungi, y = 0, yend = 10000),
             color = col_microfungi, linewidth = 0.5) +
      # median annotation
      annotate("text", x = 5e1, y = 16000,
           label = paste0("median: ", format(median_raw, big.mark = ",")),
           color = col_raw, size = 3, hjust = 0) +
      annotate("text", x = 2e1, y = 9000,
           label = paste0("median: ", format(median_microfungi, big.mark = ",")),
           color = col_microfungi, size = 3, hjust = 0) +
      theme(panel.border = element_rect(color = "black", linewidth = 1),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.position = c(0.25, 0.25),
            legend.key.size = unit(0.15, "in"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "his_depth2.pdf"), his_depth2, width = 3, height = 3)

# --- Visualize proportions of read types ---
# keep samples with analyzable reads >= 10,000
depth_filter <- depth %>% filter(depth_all >= 10000)

nrow(depth) # 47300
nrow(depth_filter) # 42366

# compute per-sample relative abundance using depth_all as denominator
depth_filter$ra_microfungi <- (depth_filter$depth_microfungi)/depth_filter$depth_all
depth_filter$ra_macrofungi <- (depth_filter$depth_macrofungi)/depth_filter$depth_all
depth_filter$ra_plant <- (depth_filter$depth_plant)/depth_filter$depth_all
depth_filter$ra_other <- (depth_filter$depth_all - depth_filter$depth_microfungi - depth_filter$depth_macrofungi - depth_filter$depth_plant)/depth_filter$depth_all

# compute overall proportion per type
ratio_microfungi <- round(sum(depth_filter$depth_microfungi) / sum(depth_filter$depth_all) * 100, 1)
ratio_macrofungi <- round(sum(depth_filter$depth_macrofungi) / sum(depth_filter$depth_all) * 100, 1)
ratio_plant <- round(sum(depth_filter$depth_plant) / sum(depth_filter$depth_all) * 100, 1)
ratio_other <- round(
  (sum(depth_filter$depth_all) -
   sum(depth_filter$depth_microfungi) -
   sum(depth_filter$depth_macrofungi) -
   sum(depth_filter$depth_plant)
  ) / sum(depth_filter$depth_all) * 100, 1)

ratio_microfungi
ratio_macrofungi
ratio_plant
ratio_other

# select the col for p
depth_filter2 <- depth_filter %>% dplyr::select(ID, PRJ, ra_microfungi, ra_macrofungi, ra_plant, ra_other)
head(depth_filter2)

depth_filter2_long <- depth_filter2 %>%
      pivot_longer(cols = c("ra_microfungi", "ra_macrofungi", "ra_plant", "ra_other"),
               names_to = "category",
               values_to = "ra")
head(depth_filter2_long)


color <- c("ra_other" = "#76877d",
           "ra_plant" = "#fff8f0",
           "ra_macrofungi" = "#92140c",
           "ra_microfungi" = "#1e1e24")
label <- c(
  ra_other = paste0("other (", ratio_other, "%)"),
  ra_plant = paste0("plant (", ratio_plant, "%)"),
  ra_macrofungi = paste0("macrofungi (", ratio_macrofungi, "%)"),
  ra_microfungi = paste0("microfungi (", ratio_microfungi, "%)")
)
# order samples by ra_microfungi (desc)
ordered_IDs <- depth_filter2_long %>%
  filter(category == "ra_microfungi") %>%
  arrange(desc(ra)) %>%
  pull(ID) %>%
  unique()

# assign sequential sample index ID2
depth_filter2_long <- depth_filter2_long %>%
  mutate(ID2 = match(ID, ordered_IDs))

# set custom stacking order
depth_filter2_long <- depth_filter2_long %>%
  mutate(category = factor(category, levels = c("ra_other",
                                                "ra_plant",
                                                "ra_macrofungi",
                                                "ra_microfungi")))
head(depth_filter2_long)

# plot
bar_reads <- ggplot(depth_filter2_long, aes(x = ID2,y = ra, fill = category)) +
             geom_bar(stat = "identity", position = "stack", width = 1, alpha = 0.85) +  # width = 1 removes gaps
             scale_fill_manual(values = color,
                               labels = label,
                               breaks = c("ra_microfungi", "ra_macrofungi", "ra_plant", "ra_other")) +
             labs(x = "Samples (n = 42,366, clean reads >= 10,000)",
                  y = "Relative abundance",
                  fill = NULL) +
             scale_x_continuous(breaks = seq(0, max(depth_filter2_long$ID2), by = 5000)) +
             scale_x_continuous(expand = c(0, 0)) + # remove axis expansion (no gap to border)
             scale_y_continuous(expand = c(0, 0)) +
             theme_bw() +
             theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.text.x = element_blank(), # hide x-axis tick labels
  axis.ticks.x = element_blank(),
  #panel.border = element_rect(color = "black", linewidth = 1),
  legend.position = c(0.2, 0.25),
  legend.key.size = unit(0.15, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_reads.png"), bar_reads, width = 6, height = 3, dpi = 900) # save PNG to avoid color issues
#ggsave(file.path(path, "bar_reads.pdf"), bar_reads, width = 6, height = 3)


# --- Diversity distribution: all species vs. microfungi ---
diversity_all <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check/div_all.csv", row.names = 1)
diversity_microfungi <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Main/data/diversity/diversity_microfungi_filterdp.csv", row.names = 1)

colnames(diversity_all) <- paste0(colnames(diversity_all), "_all")
colnames(diversity_microfungi) <- paste0(colnames(diversity_microfungi), "_microfungi")

head(diversity_all)
head(diversity_microfungi)

diversity <- merge(diversity_all, diversity_microfungi, by = 0, all = FALSE)
colnames(diversity)[1] <- "ID"

nrow(diversity_all)
nrow(diversity_microfungi)
nrow(diversity) # 37417


# Obs feature
# subset columns
diversity3 <- diversity %>% dplyr::select(ID, Observed_all, Observed_microfungi)
head(diversity3)

# compute medians
median_all <- round(median(diversity3$Observed_all), 2)
median_microfungi <- round(median(diversity3$Observed_microfungi),2)
median_all
median_microfungi

# convert to long format
diversity_long3 <- pivot_longer(
  diversity3,
  cols = c(Observed_all, Observed_microfungi),
  names_to = "category",
  values_to = "index")

diversity_long3

col <- c("Observed_all" = col_all,
         "Observed_microfungi" = col_microfungi)

label <- c("Observed_all" = "all species",
           "Observed_microfungi" = "microfungi")

his_diversity3 <- ggplot(diversity_long3, aes(x = index, fill = category)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.5, color = "black", linewidth = 0.25) +
      scale_fill_manual(values = col,
                        labels = label,
                        breaks = c("Observed_all", "Observed_microfungi")) + # fixed order
      scale_x_log10() +
      #coord_cartesian(xlim = c(1, 300)) +
      labs(x = "Observed features",
           y = "Number of samples",
           fill = NULL) +
      theme_bw() +
      # median line
      geom_segment(aes(x = median_all, xend = median_all, y = 0, yend = 4000),
             color = "#ad6a6c", linewidth = 0.5) +
      geom_segment(aes(x = median_microfungi, xend = median_microfungi, y = 0, yend = 4000),
             color = col_microfungi, linewidth = 0.5) +
      # median annotation
      annotate("text", x = 150, y = 2500,
           label = paste0("median: ", format(median_all, big.mark = ",")),
           color = "#ad6a6c", size = 3, hjust = 0) +
      annotate("text", x = 1.5, y = 2500,
           label = paste0("median: ", format(median_microfungi, big.mark = ",")),
           color = col_microfungi, size = 3, hjust = 0) +
      theme(panel.border = element_rect(color = "black", linewidth = 1),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.position = c(0.23, 0.85),
            legend.key.size = unit(0.15, "in"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "his_diversity3.pdf"), his_diversity3, width = 3, height = 3)






# --- Proportion of taxa with prevalence > threshold ---
prev_all_level_steps <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/prev_all_level_steps.csv")

col <- c("Phylum" = "#001427",
         "Class" = "#708d81",
         "Order" = "#f4d58d",
         "Family" = "#274c77",
         "Genus" = "#8d0801")

line_prev <- ggplot(prev_all_level_steps, aes(x = prevalence_threshold, y = taxa_proportion, color = level)) +
  geom_line(linewidth = 0.75, alpha = 0.8) +
  scale_color_manual(values = col,
                    breaks = c("Phylum","Class","Order","Family","Genus")) + # fixed order
  scale_y_log10()+
  labs(x = "Prevalence",
       y = "Proportion of taxa",
       color = NULL) + # no legend title
guides(color = guide_legend(nrow = 2)) +  # legend in two rows
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.55, 0.9),
        legend.key.size = unit(0.1, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "line_prev.pdf"), line_prev, width = 3, height = 3)


# --- Top 20 taxa per taxonomic level ---
prev <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/prev_all_level.csv", row.names = 1)
head(prev)

# add column and trim prefixes in taxon names
prev$taxon <- sub("^[a-z]{1,2}__", "", row.names(prev))

# append level suffix to 'Unassigned'
prev$taxon <- ifelse(prev$taxon == "Unassigned",
                   paste0("Unassigned_", prev$level),
                   prev$taxon)

# select top 20 taxa by prevalence per level
top_prev <- prev %>%
  group_by(level) %>%
  slice_max(order_by = prevalence, n = 20) %>%
  ungroup() %>%
  # set level order
  mutate(level = factor(level, levels = c("Phylum", "Class", "Order", "Family", "Genus"))) %>%
  group_by(level) %>%
  # order taxa by prevalence (desc)
  mutate(taxon = factor(taxon, levels = taxon[order(prevalence)])) %>%
  ungroup()


bar_prev <- ggplot(top_prev, aes(x = prevalence, y = taxon, fill = level)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = col) +
  geom_text(aes(x = 0.02, label = taxon),
            hjust = 0, size = 2.5) +  # show taxon labels to the right of bars
  facet_wrap(~ level, scales = "free_y", nrow = 1) +  # single-row layout
  labs(x = "Prevalence", y = "Top 20 taxa", fill = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # tilt x labels 45°
        axis.text.y = element_blank(),        # hide y-axis labels
        axis.ticks.y = element_blank(),       # hide y-axis ticks
        strip.text = element_text(face = "bold"),
        legend.position = "none" )

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_prev.pdf"), bar_prev, width = 6, height = 6)



# --- Relative abundance by taxonomic level (top 20 taxa) ---
# top 20 taxa per level
abund <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/abund_all_level.csv", row.names = 1)
head(abund)

# add column and trim prefixes in taxon names
abund$taxon <- sub("^[a-z]{1,2}__", "", row.names(abund))

# append level suffix to 'Unassigned'
abund$taxon <- ifelse(abund$taxon == "Unassigned",
                   paste0("Unassigned_", abund$level),
                   abund$taxon)

# select top 20 taxa by abundance per level
top20_abund <- abund %>%
  group_by(level) %>%
  slice_max(order_by = abundance, n = 20) %>%
  ungroup() %>%
  # set level order
  mutate(level = factor(level, levels = c("Phylum", "Class", "Order", "Family", "Genus"))) %>%
  group_by(level) %>%
  # order taxa by abundance (desc)
  mutate(taxon = factor(taxon, levels = taxon[order(abundance)])) %>%
  ungroup()


bar_abund <- ggplot(top20_abund, aes(x = abundance, y = taxon, fill = level)) +
  geom_col(alpha = 0.7) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +  # ~4 pretty x-axis ticks
  scale_fill_manual(values = col) +
  geom_text(aes(x = 0, label = taxon),
            hjust = 0, size = 2.5) +  # show taxon labels at bar start
  facet_wrap(~ level, scales = "free", nrow = 1) +  # free x and y scales per facet
  labs(x = "Relative abundance", y = "Top 20 taxa", fill = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # tilt x labels 45°
        axis.text.y = element_blank(),        # hide y-axis labels
        axis.ticks.y = element_blank(),       # hide y-axis ticks
        strip.text = element_text(face = "bold"),
        legend.position = "none" )

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_abund.pdf"), bar_abund, width = 6, height = 6)






# --- Rarefaction curve: unique taxa vs. sample count ---
raref <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/raref_all_level.csv", row.names = 1)
head(raref)
nrow(raref)

# helper: custom scientific notation formatter
sci_format_custom <- function(digits = 1) {
  function(x) {
    s <- format(x, scientific = TRUE, digits = digits)
    s <- gsub("e\\+0?(\\d+)", "e+\\1", s)
    s <- gsub("e\\-0?(\\d+)", "e-\\1", s)
    s
  }
}

line_raref <- ggplot(raref, aes(x = sample_size, y = unique_taxa, color = level)) +
  # mean line
  stat_summary(fun = mean, geom = "line", linewidth = 0.75) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 100) +  # error bars (±1 SD)
  coord_cartesian(xlim = c(0, 41000)) +
  scale_x_continuous(labels = sci_format_custom(digits = 1)) +
  scale_y_continuous(breaks = seq(0, 2500, by = 500)) +
  scale_color_manual(values = col,
                    breaks = c("Phylum","Class","Order","Family","Genus")) + # fixed order
  labs(x = "Samples",
       y = "Unique taxa",
       color = NULL) + # no legend title
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.75, 0.6))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "line_raref_sample.pdf"), line_raref, width = 3, height = 3)
