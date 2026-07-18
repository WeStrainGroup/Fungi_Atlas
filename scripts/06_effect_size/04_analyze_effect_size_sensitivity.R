# Purpose: Test balanced continental sampling, technical-factor adjustment, and Unknown-genus removal in effect-size analyses.
# Panels: Figures S7A-C.

library(tidyverse)
library(vegan)
library(broom)
library(effectsize)

# --- Construct balanced continental sample sets ---
setwd("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis")
meta <- read.csv("ID_meta_final.csv", header = T, row.names = 1)
otu_gen_dist <- readRDS("otu_gen_0.001_robait_filterdp_dist.rds")
arc_dis <- as.matrix(otu_gen_dist)

sampled_id <- function(i){
  set.seed(123 + i)
  sample_size_per_continent <- table(meta$Continent) %>% min()
  id <- meta %>%
    select(ID, Continent, Country) %>%
    group_by(Continent) %>%
    slice_sample(n = sample_size_per_continent) %>%
    ungroup() %>%
    pull(ID) %>%
    intersect(colnames(arc_dis))
  return(data.frame(ID = id, sampled_times = i))
}
sampled_id_list <- map_dfr(1:100, sampled_id)
write.csv(sampled_id_list, "sensitive_analysis/subsample_id.csv", row.names = F)


# --- Alpha-diversity effect sizes in balanced samples ---
div <- read.csv("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/diversity_microfungi_filterdp.csv", row.names = 1) %>%
  rownames_to_column("ID")
meta <- read.csv("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/ID_meta_final.csv", header = T, row.names = 1)
get_eta_result_for_subsample <- function(sampled_times){
  id_list <- read.csv("sensitive_analysis/subsample_id.csv") %>%
    filter(.data$sampled_times == .env$sampled_times) %>%
    pull(ID)

  sub_div <- div %>% filter(ID %in% id_list)
  div_meta <- sub_div %>% left_join(meta, by = "ID")

  get_eta_function <- function(var){
    formula_lm <- as.formula(paste0("Shannon ~ ", var))
    lm_model <- lm(formula_lm, data = div_meta, na.action = na.omit)

    gl <- broom::glance(lm_model)
    eta <- effectsize::eta_squared(lm_model, partial = FALSE)

    df_kruskal <- div_meta %>%
      dplyr::select(Shannon, !!sym(var)) %>%
      na.omit()

    result <- data.frame(
      diversity_index = "Shannon",
      factor = var,
      n_sample = gl$df.residual + 1,
      r_squared = gl$r.squared,
      r_squared_adj = gl$adj.r.squared,
      eta_squared = eta$Eta2[eta$Parameter == var],
      p_lm = gl$p.value,
      p_anova = anova(lm_model)$`Pr(>F)`[1],
      p_kruskal = kruskal.test(formula_lm, data = df_kruskal)$p.value,
      stringsAsFactors = FALSE
    )
    return(result)
  }
  vars_to_test <- c("DNA_extract_kit", "DNA_extract_lyticase", "DNA_extract_bead_beating",
                    "PCR_marker", "PCR_primer_pair_name", "PCR_enzyme",
                    "Seq_platform", "Seq_cycle", "Depth_microfungi_group",
                    "Sample_type", "Sex", "Age_group",
                    "Continent", "Country", "Region",
                    "PRJ")
  res <- map_dfr(vars_to_test, get_eta_function) %>% mutate(subsample = sampled_times) %>% relocate(subsample, .after = diversity_index)
  return(res)
}
final_result <- map_dfr(1:10, get_eta_result_for_subsample)
result <- final_result %>%
  group_by(subsample) %>%
  mutate(p_lm_adj = p.adjust(p_lm, method = "fdr"),
         p_anova_adj = p.adjust(p_anova, method = "fdr"),
         p_kruskal_adj = p.adjust(p_kruskal, method = "fdr")) %>%
  ungroup()
write.csv(result, "sensitive_analysis/effect_size_diversity_lm_subsample.csv", row.names = F)



# --- PERMANOVA effect sizes in balanced samples ---
get_permanova_res <- function(sampled_times){
  id_list <- read.csv("sensitive_analysis/subsample_id.csv") %>%
    filter(.data$sampled_times == .env$sampled_times) %>%
    pull(ID)
  otu_gen_dist <- readRDS("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/otu_gen_0.001_robait_filterdp_dist.rds")
  sub_otu_gen_dist <- as.matrix(otu_gen_dist)[id_list, id_list] %>% as.dist()

  meta <- read.csv("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/ID_meta_final.csv", header = T, row.names = 1)
  meta <- meta[match(id_list, meta$ID), ]

  adonis <- function(factor0){
    formula <- as.formula(paste("sub_otu_gen_dist ~", factor0))
    adonis2 <- adonis2(formula, data = meta, permutations = how(nperm = 999), na.action = na.omit, parallel = 3)
    res <- data.frame(subsample = sampled_times, factor = factor0, p.value = adonis2$`Pr(>F)`[1], R2 = adonis2$R2[1])
    return(res)
  }
  factor_list <- readLines("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/PERMANOVA/subsample/factor_list.txt")
  res_permanova_subsample <- map_dfr(factor_list, adonis) %>%
    mutate(sampled_times = sampled_times)
}
result <- map_dfr(1:10, get_permanova_res)
write.table(result, "sensitive_analysis/res_permanova_subsample.tsv", sep = "\t", row.names = F, quote = F)


# --- Sequential PERMANOVA adjusted for individual technical factors ---
get_permanova_res <- function(factor0, adj_factor0){
  otu_gen_dist <- readRDS("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/otu_gen_0.001_robait_filterdp_dist.rds")
  meta <- read.csv("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/ID_meta_final.csv", header = T, row.names = 1)
  meta <- as.data.frame(meta)
  meta <- meta[match(rownames(as.matrix(otu_gen_dist)), meta$ID), ]
  formula <- as.formula(paste("otu_gen_dist ~ ", adj_factor0, " + ", factor0))
  adonis2 <- adonis2(formula, data = meta, permutations = how(nperm = 999), na.action = na.omit, parallel = 8, by = "terms")
  saveRDS(adonis2, paste0("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/PERMANOVA/adjust/result5.0/", adj_factor0, "_", factor0, ".rds"))

  res <- data.frame(adonis2) %>% mutate(adj_factor = adj_factor0, factor = factor0)
  return(res)
}
factor_list <- c("Continent", "Country", "Region")
adj_factor_list <- c("DNA_extract_kit", "PCR_enzyme", "PCR_primer_pair_name", "Seq_platform", "PCR_marker", "Depth_microfungi_group", "DNA_extract_lyticase", "DNA_extract_bead_beating", "Seq_cycle")
param_grid <- expand_grid(
  factor0 = factor_list,
  adj_factor0 = adj_factor_list
  )
result <- pmap_dfr(param_grid, get_permanova_res)
write.table(result, "sensitive_analysis/res_permsnova_adjusttech.tsv", sep = "\t", row.names = F, quote = F)


# --- PERMANOVA sensitivity analysis after removing Unknown ---
get_permanova_res <- function(factor0){
  otu_gen_dist <- readRDS(paste0("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/unknown_check/dist_gen_dropunknown.rds"))
  meta <- read.csv("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/ID_meta_final.csv", header = T, row.names = 1)
  meta <- meta[match(rownames(as.matrix(otu_gen_dist)), meta$ID), ]
  formula <- as.formula(paste("otu_gen_dist ~", factor0))
  adonis2 <- adonis2(formula, data = meta, permutations = how(nperm = 999), na.action = na.omit, parallel = 4)

  res <- data.frame(diversity_index = "dropunknown", factor = factor0, p.value = adonis2$`Pr(>F)`[1], R2 = adonis2$R2[1])
  return(res)
}
factor_list <- readLines("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/PERMANOVA/subsample/factor_list.txt")
result <- map_dfr(factor_list, get_permanova_res)
write.table(result, "sensitive_analysis/res_permanova_dropunknown.tsv", sep = "\t", row.names = F, quote = F)
