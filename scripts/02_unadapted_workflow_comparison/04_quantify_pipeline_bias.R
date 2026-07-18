# Purpose: Quantify pipeline deviations in microbial metrics and test downstream fungus-metabolite associations.
# Panels: Figure S3B-I.

library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(ggpubr)
library(ggridges)
library(broom)
library(parallel)
library(ggnewscale)
library(ComplexUpset)
library(vegan)


cols <- c(
  "dissimilarity" = "#708d81",
  "abundance" = "grey30",
  "diversity" = "#f4d58d"
)


# - data input -
project <- "hmp"
ps_mg <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/HMP/3_phyloseq/microfungi/ASV/HMP_ps_microfungi_filterdp.rds")
ps_ua1 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted1/HMP/3_phyloseq/all/HMP_ps_all.rds")
ps_ua2 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted2/HMP/3_phyloseq/microfungi/ASV/HMP_ps_microfungi_filterdp.rds")

project <- "wep1"
ps_mg <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/WeP1/3_phyloseq/microfungi/ASV/WeP1_ps_microfungi_filterdp.rds")
ps_ua1 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted1/WeP1/3_phyloseq/all/WeP1_ps_all.rds")
ps_ua2 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted2/WeP1/3_phyloseq/microfungi/ASV/WeP1_ps_microfungi_filterdp.rds")

project <- "chgm"
ps_mg <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/CHGM/3_phyloseq/microfungi/ASV/CHGM_ps_microfungi_filterdp.rds")
ps_ua1 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted1/CHGM/3_phyloseq/all/CHGM_ps_all.rds")
ps_ua2 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted2/CHGM/3_phyloseq/microfungi/ASV/CHGM_ps_microfungi_filterdp.rds")




div_mg <- estimate_richness(ps_mg, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
div_ua1 <- estimate_richness(ps_ua1, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
div_ua2 <- estimate_richness(ps_ua2, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
head(div_mg)
head(div_ua1)
head(div_ua2)


div_all <- bind_rows(div_mg %>% rownames_to_column("sample") %>% mutate(method = "mycogap"),
                     div_ua1 %>% rownames_to_column("sample") %>% mutate(method = "unadapted1"),
                     div_ua2 %>% rownames_to_column("sample") %>% mutate(method = "unadapted2")) %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson), names_to = "metrics", values_to = "value") %>%
  mutate(cohort = project, type = "diversity")

head(div_all)



div_all_filter <- div_all %>%
  group_by(sample) %>%
  filter(n_distinct(method) == 3) %>%
  ungroup()

nrow(div_all) # 2010
nrow(div_all_filter) # 1728



ps_mg <- transform_sample_counts(ps_mg, function(x) x / sum(x))
ps_mg_gen <- aggregate_taxa(ps_mg, level = "Genus")
otu_mg_gen <- otu_table(ps_mg_gen) %>% as.data.frame()
otu_mg_gen[1:3, 1:3]

ps_ua1 <- transform_sample_counts(ps_ua1, function(x) x / sum(x))
ps_ua1_gen <- aggregate_taxa(ps_ua1, level = "Genus")
otu_ua1_gen <- otu_table(ps_ua1_gen) %>% as.data.frame()
otu_ua1_gen[1:3, 1:3]

ps_ua2 <- transform_sample_counts(ps_ua2, function(x) x / sum(x))
ps_ua2_gen <- aggregate_taxa(ps_ua2, level = "Genus")
otu_ua2_gen <- otu_table(ps_ua2_gen) %>% as.data.frame()
otu_ua2_gen[1:3, 1:3]



to_long <- function(otu, method_name) {
    otu %>%
    rownames_to_column("metrics") %>%
    pivot_longer(
      cols = -metrics,
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(method = method_name)
}

abun_all <- bind_rows(
  to_long(otu_mg_gen, "mycogap"),
  to_long(otu_ua1_gen, "unadapted1"),
  to_long(otu_ua2_gen, "unadapted2")) %>%
mutate(cohort = project, type = "abundance")

head(abun_all)
nrow(abun_all)



metrics_method_count <- abun_all %>%
  group_by(metrics) %>%
  summarise(n_method = n_distinct(method)) %>%
  pull(n_method)

table(metrics_method_count)


abun_all_filter <- abun_all %>%

  group_by(sample) %>%
  filter(n_distinct(method) == 3) %>%
  ungroup() %>%


  group_by(sample, metrics) %>%
  filter(sum(value) > 0) %>%
  ungroup()

nrow(abun_all) # 485080
nrow(abun_all_filter) # 4761





dist_mg <- vegdist(otu_mg_gen %>% t(), method = "bray")
dist_ua1 <- vegdist(otu_ua1_gen %>% t(), method = "bray")
dist_ua2 <- vegdist(otu_ua2_gen %>% t(), method = "bray")
str(dist_mg)


dist_to_long <- function(dist_obj, method_name) {
    as.matrix(dist_obj) %>%
    as.data.frame() %>%
    rownames_to_column("sample1") %>%
    pivot_longer(
      cols = -sample1,
      names_to = "sample2",
      values_to = "value"
    ) %>%
    filter(sample1 < sample2) %>%
    mutate(sample = paste0(sample1, "_", sample2),
           method = method_name,
           metrics = "Bray-Curtis") %>%
    select(sample, method, metrics, value)
}



dist_all <- bind_rows(
  dist_to_long(dist_mg, "mycogap"),
  dist_to_long(dist_ua1, "unadapted1"),
  dist_to_long(dist_ua2, "unadapted2")
) %>%
  mutate(cohort = project, type = "dissimilarity")

head(dist_all)



dist_all_filter <- dist_all %>%
  group_by(sample) %>%
  filter(n_distinct(method) == 3) %>%
  ungroup()

nrow(dist_all) # 77427
nrow(dist_all_filter) # 55008





mm_all <- rbind(div_all_filter, abun_all_filter, dist_all_filter)
head(mm_all)
tail(mm_all)
nrow(mm_all)
table(is.na(mm_all$value))


mm_all <- mm_all %>%
  group_by(metrics) %>%
  mutate(value_scale = (value - min(value, na.rm = TRUE)) / (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
  mutate(value_scale = ifelse(is.na(value_scale), 1, value_scale)) %>%
  ungroup() %>%
  as.data.frame()

head(mm_all)
nrow(mm_all)
table(is.na(mm_all$value_scale))



mm_all_wide_noscale <- mm_all %>%
    select(sample, method, metrics, type, value) %>%
    pivot_wider(names_from = method, values_from = value, values_fill = 0) %>%
    mutate(type = factor(type, levels = c("dissimilarity", "abundance", "diversity"))) %>%
    arrange(sample) %>%
    as.data.frame()

mm_all_wide <- mm_all %>%
    select(sample, method, metrics, type, value_scale) %>%
    pivot_wider(names_from = method, values_from = value_scale, values_fill = 0) %>%
    mutate(type = factor(type, levels = c("dissimilarity", "abundance", "diversity"))) %>%
    arrange(sample) %>%
    as.data.frame()

head(mm_all_wide)
tail(mm_all_wide)





scatter <- ggplot(mm_all_wide, aes(x = mycogap, y = unadapted1, color = type)) +
  geom_point(alpha = 0.5, shape = 16, size = 2) +
  scale_color_manual(values = cols) +
  geom_abline(
    slope = 1, intercept = 0,
    color = "black",
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  labs(
    x = "Scaled MMs of MycoGAP",
    y = "Scaled MMs of unadapted pipeline"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
  legend.position = "none")

path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
ggsave(file.path(path, paste0("scatter_", project, ".pdf")), scatter, width = 3, height = 3)






mm_bias <- mm_all_wide %>%
  group_by(type) %>%
  mutate(

    min_nonzero = min(c(mycogap[mycogap > 0],
                        unadapted1[unadapted1 > 0],
                        unadapted2[unadapted2 > 0])),


    mycogap_adj = mycogap + (min_nonzero / 2),
    unadapted1_adj = unadapted1 + (min_nonzero / 2),
    unadapted2_adj = unadapted2 + (min_nonzero / 2),

    # bias
    unadapted1_bias = abs(mycogap_adj - unadapted1_adj) / mycogap_adj,
    unadapted2_bias = abs(mycogap_adj - unadapted2_adj) / mycogap_adj
  ) %>%
  ungroup() %>%
  select(sample, metrics, type, unadapted1_bias, unadapted2_bias) %>%
  as.data.frame()

head(mm_bias)
tail(mm_bias)



mm_bias_plot <- mm_bias %>%
  mutate(unadapted1_bias_plot = pmin(pmax(unadapted1_bias, 1e-4), 1e4),
         unadapted2_bias_plot = pmin(pmax(unadapted2_bias, 1e-4), 1e4))

min(mm_bias_plot$unadapted1_bias_plot)
max(mm_bias_plot$unadapted1_bias_plot)


prop_df <- mm_bias_plot %>%
  group_by(type) %>%
  summarise(prop_gt1 = mean(unadapted1_bias > 1, na.rm = TRUE)) %>%
  mutate(label = paste0("> ", "1: ", round(prop_gt1 * 100, 1), "%"))
prop_df


his <- ggplot(mm_bias_plot, aes(x = unadapted1_bias_plot, y = type, fill = type)) +
  geom_density_ridges(
    alpha = 0.5,
    scale = 1,
    color = "black",
    rel_min_height = 0,
    bandwidth = 0.2) +
  scale_x_log10() +
  scale_fill_manual(values = cols) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5) +
  geom_text(
    data = prop_df,
    aes(x = 10, y = type, label = label),
    hjust = 0,
    nudge_y = 0.25,
    size = 3) +
  geom_text(
    data = prop_df,
    aes(x = 0.000025, y = type, label = type),
    hjust = 0,
    nudge_y = 0.75,
    size = 3) +
  labs(
    x = "Relative bias of MMs",
    y = "Density",
    fill = NULL,
    color = NULL) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.key.size = unit(0.15, "in"),
            legend.position = "none")


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
ggsave(file.path(path, paste0("his_", project, ".pdf")), his, width = 3, height = 3)





mm_bias_stats <- mm_bias %>%
  group_by(type) %>%
  summarise(
    mean_unadapted1 = mean(unadapted1_bias, na.rm = TRUE),
    sd_unadapted1   = sd(unadapted1_bias, na.rm = TRUE),
    mean_unadapted2 = mean(unadapted2_bias, na.rm = TRUE),
    sd_unadapted2   = sd(unadapted2_bias, na.rm = TRUE),
  ) %>%
  mutate(mean_ratio = mean_unadapted1 / mean_unadapted2,
        improvement_pct = (mean_unadapted1 - mean_unadapted2) / mean_unadapted1 * 100) %>%
  as.data.frame()

mm_bias_stats


mm_bias_long <- mm_bias %>%
  pivot_longer(cols = c(unadapted1_bias, unadapted2_bias),
               names_to = "bias_type",
               values_to = "bias_value") %>%
  group_by(type, bias_type) %>%
  mutate(min_nonzero = min(bias_value[bias_value > 0], na.rm = TRUE),
         bias_value_plot = ifelse(bias_value == 0, min_nonzero / 2, bias_value)) %>%
  ungroup() %>%
  select(-min_nonzero) %>%
  mutate(type = factor(type, levels = c("diversity", "abundance", "dissimilarity")))


head(mm_bias_long)
mean()



box <- ggplot(mm_bias_long, aes(x = bias_type, y = bias_value_plot, fill = bias_type)) +
  geom_boxplot(alpha = 0.7, color = "black", outlier.shape = 16) +
  scale_x_discrete(labels = c("unadapted1_bias" = "Unadapted", "unadapted2_bias" = "Microfungi only")) +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.11))) +
  facet_wrap(~type, scales = "free_y") +
  theme_bw() +
  labs(y = "Relative bias of MMs") +
  scale_fill_manual(values = c("unadapted1_bias" = "grey85", "unadapted2_bias" = "#606c38")) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = list(c("unadapted1_bias","unadapted2_bias")),
    label = "p.format",
    size = 3,
    tip.length = 0,
    bracket.size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
ggsave(file.path(path, paste0("box_", project, ".pdf")), box, width = 4, height = 3)





path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
write.csv(mm_all_wide_noscale, file.path(path, paste0("mm_data_", project, ".csv")))
write.csv(mm_bias, file.path(path, paste0("mm_bias_", project, ".csv")))
write.csv(mm_bias_stats, file.path(path, paste0("mm_bias_stats_", project, ".csv")))
















codebook <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Omics/WeP1_fecal_metabolite_codebook.csv", row.names = 1)
metabo <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Omics/WeP1_fecal_metabolite.csv", row.names = 1)
meta <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Omics/WeP1_pheno_core.csv", row.names = 1)

# metadata
meta <- meta %>%
    select(age, sex) %>%
    mutate(sex = factor(sex)) %>%
    rownames_to_column("sample") %>%
    mutate(sample = paste0("V1", sample))

head(meta)


metabo <- metabo %>%
  rownames_to_column("sample") %>%
  mutate(sample = str_replace(
    sample,
    "^(WPN[0-9]+)(V[0-9]+).$",
    "\\2\\1"
  )) %>%
  filter(str_detect(sample, "^V1")) %>%
  column_to_rownames("sample") %>%
  select(where(~ !all(. == 0)))

metabo[1:3, 1:3]
dim(metabo) # 170 * 970


codebook_filter <- codebook %>% filter(Level != 3)
nrow(codebook) # 1815
nrow(codebook_filter) # 1258




min(apply(metabo, 2, min))

prev <- colMeans(metabo != 9)
prev_high <- names(prev[prev >= 0.5])
length(prev_high) # 960



keep_cols <- intersect(prev_high, rownames(codebook_filter))
length(keep_cols)

metabo_filter <- metabo %>%
    select(all_of(keep_cols))
dim(metabo_filter) # 170 * 960



metabo_filter_trans <- metabo_filter %>%
  mutate(across(where(is.numeric), log)) %>%
  mutate(across(where(is.numeric), ~ as.numeric(scale(.))))

# check
metabo_filter_trans[1:3, 1:3]
apply(metabo_filter, 2, mean)
apply(metabo_filter, 2, sd)

apply(metabo_filter_trans, 2, mean)
apply(metabo_filter_trans, 2, sd)




myco <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis/mm_data_wep1.csv", row.names = 1)

myco <- myco %>%
    filter(type != "dissimilarity") %>%
    filter(str_detect(sample, "^V1")) %>%
    filter(!(mycogap == 0 & unadapted1 == 0 & unadapted2 == 0))

head(myco)


ps <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/WeP1/3_phyloseq/microfungi/ASV/WeP1_ps_microfungi_filterdp.rds")
ps <- subset_samples(ps, grepl("^V1", sample_names(ps)))
ps_gen <- aggregate_taxa(ps_mg, level = "Genus")

prev <- as.data.frame(prevalence(ps_mg_gen))
colnames(prev) <- "prev"
taxon_to_keep <- prev %>%
  filter(prev > 0.1) %>%
  rownames() %>%
  .[!grepl("Unknown", .)]
length(taxon_to_keep)


metrics_to_keep <- c(taxon_to_keep, "Observed", "Shannon", "Simpson")
metrics_to_keep


myco_long <- myco %>%
  filter(metrics %in% metrics_to_keep) %>%
  pivot_longer(cols = c(mycogap, unadapted1, unadapted2),
               names_to = "pipeline", values_to = "value") %>%
  mutate(var = paste(metrics, pipeline, sep = "_"))

head(myco_long)


myco_wide <- myco_long %>%
  select(sample, var, value) %>%
  pivot_wider(names_from = var, values_from = value, values_fill = 0) %>%
  column_to_rownames("sample")

myco_wide[1:5, 1:5]


var_type <- myco_long %>% distinct(var, type)
abund_cols <- var_type$var[var_type$type == "abundance"]
div_cols  <- var_type$var[var_type$type == "diversity"]


min_nonzero <- mycotab %>%
  filter(type == "abundance") %>%
  select(mycogap, unadapted1, unadapted2) %>%
  unlist() %>%
  (\(x) min(x[x > 0], na.rm = TRUE))()
min_nonzero



myco_trans <- myco_wide %>%
  mutate(
    across(all_of(abund_cols), ~ {
      x <- .
      x[x == 0] <- min_nonzero / 2
      log(x)
    }),

    across(all_of(div_cols), log)
  ) %>%
   mutate(across(where(is.numeric), ~ as.numeric(scale(.))))

myco_trans[1:3, 1:3]



metabo_filter_trans[1:3, 1:3]
myco_trans[1:3, 1:3]
head(meta)

df_full <- meta %>%
  inner_join(myco_trans %>% rownames_to_column("sample"), by = "sample") %>%
  inner_join(metabo_filter_trans %>% rownames_to_column("sample"), by = "sample")  %>%
  column_to_rownames("sample")

nrow(df_full) # 148
df_full[1:4,1:4]


fungi <- c(div_cols, abund_cols)
metabo <- colnames(metabo_filter_trans)
comb <- expand.grid(fungi = fungi, metab = metabo, stringsAsFactors = FALSE)
comb



res_list <- mclapply(1:nrow(comb), function(i) {

  fungi_col <- comb$fungi[i]
  met_col <- comb$metab[i]

  cat(fungi_col, "~", met_col, "\n")

  fit <- lm(df_full[[met_col]] ~ df_full[[fungi_col]] + df_full$age + df_full$sex)

  coef_mat <- coef(summary(fit))
  gl <- broom::glance(fit)

  data.frame(
    fungi = fungi_col,
    metab = met_col,
    n_sample = gl$df.residual + 1,
    intercept = coef_mat[1,1],
    slope = coef_mat[2,1],
    r_squared = gl$r.squared,
    r_squared_adj = gl$adj.r.squared,
    p_lm = gl$p.value,
    p_fungi = coef_mat[2,4]
  )

}, mc.cores = 8)



lm_result <- do.call(rbind, res_list)
lm_result <- lm_result %>%
  separate(fungi, into = c("metrics", "method"), sep = "_(?=[^_]+$)") %>%
  group_by(method, metrics) %>%
  mutate(p_lm_adj = p.adjust(p_lm, method = "BH"),
         p_fungi_adj = p.adjust(p_fungi, method = "BH")) %>%
  ungroup() %>%
  arrange(p_fungi_adj) %>%
  as.data.frame()

head(lm_result)
tail(lm_result)
nrow(lm_result)


lm_result_filter <- lm_result %>% filter(p_fungi_adj < 0.05)

head(lm_result_filter)
nrow(lm_result_filter) # 28
length(unique(lm_result_filter$metrics)) # 11
length(unique(lm_result_filter$metab)) # 17



lm_result_plot <- lm_result %>%
        filter(metrics %in% unique(lm_result_filter$metrics)) %>%
        filter(metab %in% unique(lm_result_filter$metab)) %>%
        mutate(metrics = gsub("^g__", "", metrics)) %>%
        mutate(metrics_method = paste0(metrics, "_", method)) %>%
        mutate(beta_plot = ifelse(p_fungi_adj < 0.05, slope, 0)) %>%
        select(metrics, metab, method, beta_plot, metrics_method, p_fungi_adj) %>%
        mutate(method = factor(method))


nrow(lm_result_plot) # 561
head(lm_result_plot)
tail(lm_result_plot)




anno_data <- lm_result_plot %>% filter(beta_plot == 0)


heat <- ggplot(lm_result_plot, aes(y = metrics_method, x = metab)) +

  geom_tile(aes(fill = beta_plot)) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "beta"
  ) +


  ggnewscale::new_scale_fill() +


  geom_tile(data = anno_data, aes(fill = method), color = NA, alpha = 0.2) +
  scale_fill_manual(
    name = "pipeline",
    values = c(mycogap = "grey20", unadapted1 = "grey85", unadapted2 = "#606c38"),
    labels = c(mycogap = "MycoGAP", unadapted1 = "Unadapted", unadapted2 = "Microfungi only"),
    guide = guide_legend(override.aes = list(color = NA))) +
  facet_grid(metrics ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Fecal metabolites", y = "Microbial metrics") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    legend.key.size = unit(0.15, "in"),
    panel.spacing.y = unit(3, "pt"),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
ggsave(file.path(path, "heat_metab_fungi.pdf"), heat, width = 7, height = 3, bg = "transparent")





head(lm_result_plot)


mycogap_set <- lm_result_plot %>%
  filter(p_fungi_adj < 0.05) %>%
  filter(method == "mycogap") %>%
  mutate(pair = paste0(metrics, "_", metab)) %>%
  pull(pair)

unadapted1_set <- lm_result_plot %>%
  filter(p_fungi_adj < 0.05) %>%
  filter(method == "unadapted1") %>%
  mutate(pair = paste0(metrics, "_", metab)) %>%
  pull(pair)

unadapted2_set <- lm_result_plot %>%
  filter(p_fungi_adj < 0.05) %>%
  filter(method == "unadapted2") %>%
  mutate(pair = paste0(metrics, "_", metab)) %>%
  pull(pair)



all_pairs <- unique(c(mycogap_set, unadapted1_set, unadapted2_set))

upset_df <- data.frame(
  pair = all_pairs,
  MycoGAP = all_pairs %in% mycogap_set,
  Unadapted = all_pairs %in% unadapted1_set,
  Microfungi_only = all_pairs %in% unadapted2_set
)


upset <- upset(
  upset_df,
  intersect = c("MycoGAP", "Unadapted", "Microfungi_only"),
  name = "",
  width_ratio = 0.5,
)

path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
ggsave(file.path(path, "upset_metab_fungi.pdf"), upset, width = 4, height = 4, bg = "transparent")








dist_metabo <- vegdist(metabo_filter_trans, method = "euclidean")
dist_metabo <- as.matrix(dist_metabo)
str(dist_metabo)


ps_mg <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/WeP1/3_phyloseq/microfungi/ASV/WeP1_ps_microfungi_filterdp.rds")
ps_ua1 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted1/WeP1/3_phyloseq/all/WeP1_ps_all.rds")
ps_ua2 <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted2/WeP1/3_phyloseq/microfungi/ASV/WeP1_ps_microfungi_filterdp.rds")

ps_to_dist <- function(ps, level = "Genus", method = "bray") {
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  ps_agg <- aggregate_taxa(ps_rel, level = level)
  otu <- as.data.frame(otu_table(ps_agg))
  dist <- vegdist(otu %>% t(), method = method)
  as.matrix(dist)
}

dist_mg <- ps_to_dist(ps_mg)
dist_ua1 <- ps_to_dist(ps_ua1)
dist_ua2 <- ps_to_dist(ps_ua2)



common_samples <- Reduce(
  intersect,
  list(
    rownames(dist_metabo),
    rownames(dist_mg),
    rownames(dist_ua1),
    rownames(dist_ua2)
  )
)

common_samples

dist_metabo_filter <- dist_metabo[common_samples, common_samples]
dist_mg_filter <- dist_mg[common_samples, common_samples]
dist_ua1_filter <- dist_ua1[common_samples, common_samples]
dist_ua2_filter <- dist_ua2[common_samples, common_samples]

dist_metabo_filter[1:3, 1:3]
dist_mg_filter[1:3, 1:3]


# mantel
run_mantel <- function(d1, d2, label, method = "spearman", permutations = 999) {

  res <- mantel(d1, d2, method = method, permutations = permutations)

  tibble(
    method = label,
    n_samples = nrow(d1),
    mantel_r = unname(res$statistic),
    mantel_p = res$signif,
    mantel_method = method,
    permutations = permutations
  )
}

set.seed(6311)
mantel_res <- bind_rows(
  run_mantel(dist_metabo_filter, dist_mg_filter, label = "mycogap"),
  run_mantel(dist_metabo_filter, dist_ua1_filter, label = "unadapted1"),
  run_mantel(dist_metabo_filter, dist_ua2_filter, label = "unadapted2")
)

mantel_res


# procrustes
run_procrustes <- function(d1, d2, label, permutations = 999, k = NULL) {


  if (is.null(k)) k <- nrow(d1) - 1


  pcoa1 <- cmdscale(d1, k = k, eig = TRUE, add = TRUE)
  pcoa2 <- cmdscale(d2, k = k, eig = TRUE, add = TRUE)

  res <- protest(X = pcoa1$points, Y = pcoa2$points, permutations = permutations, symmetric = TRUE)

  tibble(
    method = label,
    n_samples = nrow(d1),
    procrustes_ss = res$ss,
    procrustes_cor = res$t0,
    procrustes_p = res$signif,
    permutations = permutations,
  )
}

set.seed(6311)
procrustes_res <- bind_rows(
  run_procrustes(dist_metabo_filter, dist_mg_filter, label = "mycogap"),
  run_procrustes(dist_metabo_filter, dist_ua1_filter, label = "unadapted1"),
  run_procrustes(dist_metabo_filter, dist_ua2_filter, label = "unadapted2")
)

procrustes_res


procrustes_res$method <- factor(procrustes_res$method,levels = c("unadapted2", "unadapted1", "mycogap"))

bar <- ggplot(procrustes_res, aes(x = procrustes_cor, y = method, fill = method)) +
  geom_col(width = 0.5, alpha = 0.7, color = "black", linewidth = 0.5) +
  labs(x = "Procrustes correlation", y = "Pipeline") +
  xlim(0, 1) +
  geom_text(aes(y = method, label = sprintf("p = %.3f", procrustes_p)), x = 0.05, hjust = 0, size = 3) +
  scale_y_discrete(labels = c(mycogap = "MycoGAP", unadapted1 = "Unadapted", unadapted2 = "Microfungi only")) +
  scale_fill_manual(name = "pipeline",
                    values = c(mycogap = "grey20", unadapted1 = "grey85", unadapted2 = "#606c38")) +
  theme_bw() +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        legend.position = "none")


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
ggsave(file.path(path, "bar_metab_fungi.pdf"), bar, width = 3, height = 3, bg = "transparent")




path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Analysis"
write.csv(lm_result_filter, file.path(path, "metab_fungi_lm.csv"))
write.csv(procrustes_res, file.path(path, "metab_fungi_procrustes.csv"))
