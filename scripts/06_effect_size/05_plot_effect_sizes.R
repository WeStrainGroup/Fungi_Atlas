# Purpose: Plot diversity, composition, and genus-level effect sizes and primer/platform differential abundance.
# Panels: Figure 3A-E.

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

# --- Effect sizes of different factors on diversity ---
div_effect <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/effect_size/effect_size_diversity_lm.csv")

# keep the target index only
div_effect <- div_effect %>% filter(div_index == "Shannon")
head(div_effect)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/effect_size/effect_size_meta.csv")
head(meta)

# merge
div_effect_meta <- merge(div_effect, meta, by = "factor", all.x = T)
head(div_effect_meta)
unique(div_effect_meta$type)

# keep factors included for analysis
div_effect_meta <- div_effect_meta %>% filter(to_analysis == 1)
head(div_effect_meta)

# add significance labels
div_effect_meta <- div_effect_meta %>%
  mutate(sig_label = ifelse(p_anova_adj > 0.01, "ns", NA))

# order factor names by eta_squared (ascending)
div_effect_meta <- div_effect_meta %>%
  mutate(name = factor(name, levels = name[order(eta_squared, decreasing = F)]))

# visualization
type_colors <- c("Demography" = "#001427",
                 "Geography" = "#708d81",
                 "Sequencing" = "#f4d58d",
                 "DNA_extract" = "#274c77",
                 "PCR" = "#8d0801",
                 "Sample_type" = "#999999")


# plot all factors together
bar_effect <- ggplot(div_effect_meta, aes(x = name, y = eta_squared, fill = type)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  coord_flip() +  # horizontal bars for readability
  #scale_fill_paletteer_d("ggsci::signature_substitutions_cosmic") +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       #y = expression("Effect size (partial "* eta^2 *")"),
       y = expression("Effect size (" * eta^2 * ")"),
       fill = NULL) +
geom_text(
    aes(label = sig_label),
    vjust = 0.5, hjust = -0.1,
    size = 3.5,
    na.rm = TRUE  # skip NA values
  ) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid.major.y = element_blank(),
        plot.margin = unit(c(5.5, 10, 5.5, 5.5), "pt"),  # margins: top, right, bottom, left
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.65, 0.25),
        legend.key.size = unit(0.15, "in"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_effect_size_div_lm.pdf"), bar_effect, width = 3, height = 3)





# --- Effect sizes of different factors on distance matrix ---
dis_effect <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/effect_size/res_PERMANOVA.csv", row.names = 1)
dis_effect <- dis_effect %>% rownames_to_column("factor")
head(dis_effect)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/effect_size/effect_size_meta.csv")
head(meta)

# merge
dis_effect_meta <- merge(dis_effect, meta, by = "factor", all.x = T)
head(dis_effect_meta)
unique(dis_effect_meta$type)

# keep factors included for analysis
dis_effect_meta <- dis_effect_meta %>% filter(to_analysis == 1)
head(dis_effect_meta)

# add significance labels
dis_effect_meta <- dis_effect_meta %>%
  mutate(sig_label = ifelse(p_value > 0.001, "ns", NA))

# order factor names by R2 (ascending)
dis_effect_meta <- dis_effect_meta %>%
  mutate(name = factor(name, levels = name[order(R2, decreasing = F)]))

# visualization
type_colors <- c("Demography" = "#001427",
                 "Geography" = "#708d81",
                 "Sequencing" = "#f4d58d",
                 "DNA_extract" = "#274c77",
                 "PCR" = "#8d0801",
                 "Sample_type" = "#999999")


# plot all factors together
bar_effect <- ggplot(dis_effect_meta, aes(x = name, y = R2, fill = type)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  coord_flip() +  # horizontal bars for readability
  #scale_fill_paletteer_d("ggsci::signature_substitutions_cosmic") +
  scale_fill_manual(values = type_colors) +
  labs(x = NULL,
       y = expression("Effect size (" * R^2 * ")"),
       fill = NULL) +
geom_text(
    aes(label = sig_label),
    vjust = 0.5, hjust = -0.1,
    size = 3.5,
    na.rm = TRUE  # skip NA values
  ) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        panel.grid.major.y = element_blank(),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.65, 0.25),
        legend.key.size = unit(0.15, "in"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_effect_size_dis.pdf"), bar_effect, width = 3, height = 3)






# --- Effects of factors on abundances of the top 25 prevalent genera ---
gen_effect <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/effect_size/effect_size_genus_lm.csv")
head(gen_effect)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/effect_size/effect_size_meta.csv")
head(meta)

# compute list of top prevalent genera
prev <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/prev_all_level.csv", row.names = 1)
head(prev)

top_genus <- prev %>%
   filter(level == "Genus", # keep Genus level
         rownames(prev) != "g__Unassigned") %>% # remove g__Unassigned
  arrange(desc(prevalence)) %>%          # order by prevalence (desc)
  slice_head(n = 25)                     # take top 25

top_genus <- rownames(top_genus)
top_genus

# merge
gen_effect_meta <- merge(gen_effect, meta, by = "factor", all.x = T)
head(gen_effect_meta)

# keep selected factors
gen_effect_meta <- gen_effect_meta %>% filter(name %in% c("Sample_type", "Region", "Age_group", "Sex", "DNA_extract", "PCR_primer", "PCR_enzyme", "Seq_platform", "Seq_cycle", "Seq_depth"))
head(gen_effect_meta)
nrow(gen_effect_meta) # 1980

# keep top-prevalence taxa
gen_effect_meta <- gen_effect_meta %>% filter(genus %in% top_genus)
head(gen_effect_meta)
nrow(gen_effect_meta) #250

# keep rows with adjusted p <= 0.01
gen_effect_meta <- gen_effect_meta %>% filter(p_anova_adj <= 0.01)
nrow(gen_effect_meta) # 231

# trim g__ prefix from genus names
gen_effect_meta$genus <- gsub("^g__", "", gen_effect_meta$genus)
head(gen_effect_meta)


# find dominant factor type per genus (max summed eta_squared)
genus_type_label <- gen_effect_meta %>%
  group_by(genus, type) %>%
  summarise(sum_eta = sum(eta_squared, na.rm = TRUE), .groups = "drop") %>%  # sum eta_squared by genus and type
  group_by(genus) %>%
  slice_max(order_by = sum_eta, n = 1, with_ties = FALSE) %>%   # keep type with max sum_eta per genus
  ungroup() %>%
  dplyr::select(genus, dominant_type = type, sum_eta)

head(genus_type_label)

# merge labels back to data
gen_effect_meta_labeled <- gen_effect_meta %>%
  left_join(genus_type_label, by = "genus")

head(gen_effect_meta_labeled)
nrow(gen_effect_meta_labeled)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/effect_size"
write.csv(gen_effect_meta_labeled, file.path(path, "effect_size_genus25_lm_plot.csv"))

# order genera by dominant_type groups
# Step 1: compute summed eta_squared for each genus × dominant_type
genus_order_df <- gen_effect_meta_labeled %>%
  group_by(genus, dominant_type) %>%
  summarise(eta2_sum = sum(eta_squared, na.rm = TRUE), .groups = "drop") %>%
  arrange(dominant_type, desc(eta2_sum))

# Step 2: set factor order
gen_effect_meta_labeled$genus <- factor(
  gen_effect_meta_labeled$genus,
  levels = genus_order_df$genus
)

# order stacked bars by dominant_type first
# Step 1: aggregate to genus × type
stacked_df <- gen_effect_meta_labeled %>%
  group_by(genus, type, dominant_type) %>%
  summarise(eta2_sum = sum(eta_squared), .groups = "drop")

# Step 2: within each genus, place dominant_type first, others by eta2_sum desc
stacked_df <- stacked_df %>%
  group_by(genus) %>%
  mutate(
    dominant_flag = ifelse(type == dominant_type, 1, 2)
  ) %>%
  arrange(dominant_flag, desc(eta2_sum), .by_group = TRUE) %>%
  mutate(
    ymin = c(0, cumsum(eta2_sum)[-n()]),
    ymax = cumsum(eta2_sum)
  ) %>%
  ungroup()

type_colors <- c("Demography" = "#001427",
                 "Geography" = "#708d81",
                 "Sequencing" = "#f4d58d",
                 "DNA_extract" = "#274c77",
                 "PCR" = "#8d0801",
                 "Sample_type" = "#999999")

bar_gen_effect <- ggplot(stacked_df, aes(x = genus, ymin = ymin, ymax = ymax, fill = type)) +
  geom_rect(aes(xmin = as.numeric(genus) - 0.4,
                xmax = as.numeric(genus) + 0.4),
            color = NA) +
  scale_fill_manual(values = type_colors) +
  labs(
    x = "Top 25 most prevalent genera",
    #y = expression("Effect size (partial "* eta^2 *")"),
    y = expression("Effect size (" * eta^2 * ")"),
    fill = NULL
  ) +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = c(0.15, 0.65),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.key.size = unit(0.15, "in")
  )

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "bar_effect_size_gen25_lm.pdf"), bar_gen_effect, width = 6, height = 3)



# Pie chart for dominant_type
# sum eta by dominant_type
pie_data <- genus_type_label %>%
  group_by(dominant_type) %>%
  summarise(sum_eta = sum(sum_eta)) %>%
  ungroup()

# compute percentages and label positions
pie_data <- pie_data %>%
  mutate(
    fraction = sum_eta / sum(sum_eta),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    label_pos = (ymax + ymin) / 2,
    label = paste0(dominant_type, "\n", scales::percent(fraction)),
    percent_label = scales::percent(fraction)
  )


pie <- ggplot(pie_data, aes(x = 2, y = sum_eta, fill = dominant_type)) +  # x=2 controls ring radius
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +  # control inner/outer radii; larger range → thicker ring
  theme_void() +
  geom_text(aes(label = percent_label), position = position_stack(vjust = 0.5), size = 4.5) +
  scale_fill_manual(values = type_colors) +
  theme(legend.position = "none")



path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "pie_effect_size_gen25_lm.pdf"), pie, width = 2, height = 2, bg = "transparent")







# --- Volcano plot: taxa enriched by different primers ---
lmm_primer <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/differential_analysis/lmm_PCR_primer.csv", row.names =1)
lmm_primer <- lmm_primer %>%
  mutate(
    taxon = str_remove(taxon, "^g__"),  # remove leading g__
    neg_log10_p = -log10(p_adj),
    sig = case_when(
      p_adj <=0.05 & estimate >0 ~ "ITS3-ITS4",
      p_adj <=0.05 & estimate <0 ~ "ITS1F-ITS2",
      TRUE ~ "Not significant"
    )
  )


p <- ggplot(lmm_primer, aes(x = estimate, y= neg_log10_p, color = sig)) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +  # FDR=0.05 threshold
    geom_point(alpha=0.75, size=2, shape = 16) +
  xlim(-2.5, 2.5) +
  ylim(0, 2.5) +
  scale_color_manual(
    values = c(
      "ITS1F-ITS2" = "#457b9d",
      "ITS3-ITS4" = "#e63946",
      "Not significant" = "grey"),
    breaks = c("ITS1F-ITS2", "ITS3-ITS4")) + # show only these two
  ggrepel::geom_text_repel(
    data = subset(lmm_primer, p_adj<=0.05),
    aes(label=taxon),
    size = 3,
    max.overlaps = 20,
    force = 1,
    color = "black",
    segment.color = "#00000040") +
  labs(
    x = "Estimated CLR abundance difference",
    y = "-log10(FDR)",
    color = " Enriched in:") +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        #panel.grid= element_blank(),
        legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(0.15, "in"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "volcano_primer.pdf"), p, width = 3, height = 3)



# --- Volcano plot: taxa differences across sequencing platforms ---
lmm_platform <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/differential_analysis/lmm_Seq_platform.csv", row.names =1)
lmm_platform <- lmm_platform %>%
  mutate(
    taxon = str_remove(taxon, "^g__"),  # remove leading g__
    neg_log10_p = -log10(p_adj),
    sig = case_when(
      p_adj <=0.05 & estimate >= 0.5 ~ "NovaSeq",
      p_adj <=0.05 & estimate <= -0.5 ~ "MiSeq",
      TRUE ~ "Not significant"
    )
  )


p <- ggplot(lmm_platform, aes(x = estimate, y= neg_log10_p, color = sig)) +
  geom_vline(xintercept = 0.5, linetype="dashed") +
  geom_vline(xintercept = -0.5, linetype="dashed") +

  geom_hline(yintercept = -log10(0.05), linetype="dashed") +  # FDR=0.05 threshold
    geom_point(alpha=0.75, size=2, shape = 16) +
  xlim(-3.5, 3.5) +
  #ylim(0, 2) +
  scale_color_manual(
    values = c(
      "MiSeq" = "#457b9d",
      "NovaSeq" = "#e63946",
      "Not significant" = "grey"),
    breaks = c("MiSeq", "NovaSeq")) + # show only these two
  ggrepel::geom_text_repel(
    data = subset(lmm_platform, p_adj<=0.05 & abs(estimate) >= 0.5),
    aes(label=taxon),
    size = 3,
    na.rm = FALSE,
    force = 10,
    color = "black",
    segment.color = "#00000040") +
  labs(
    x = "Estimated CLR abundance difference",
    y = "-log10(FDR)",
    color = " Enriched in:") +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        #panel.grid= element_blank(),
        legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.background = element_rect(fill = "#FFFFFF1A", color = NA),
        #legend.position = "bottom",
        legend.key.size = unit(0.15, "in"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/figure"
ggsave(file.path(path, "volcano_platform.pdf"), p, width = 3, height = 3)
