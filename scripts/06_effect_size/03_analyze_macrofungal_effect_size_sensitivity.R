# Purpose: Re-estimate diversity and composition effect sizes under alternative macrofungal filtering strategies.
# Panels: Figure S7A-B.

library(tidyverse)
library(vegan)
library(broom)
library(effectsize)

# --- PERMANOVA sensitivity analysis for macrofungal filtering ---
get_permanova_res <- function(shannon, factor0){
  otu_gen_dist <- readRDS(paste0("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/macrofungi_check/dist_micro_", shannon, "_gen.rds"))
  meta <- read.csv("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/ID_meta_final.csv", header = T, row.names = 1)
  meta <- meta[match(rownames(as.matrix(otu_gen_dist)), meta$ID), ]

  formula <- as.formula(paste("otu_gen_dist ~", factor0))
  adonis2 <- adonis2(formula, data = meta, permutations = how(nperm = 999), na.action = na.omit, parallel = 4)
  res <- data.frame(diversity_index = shannon, factor = factor0, p.value = adonis2$`Pr(>F)`[1], R2 = adonis2$R2[1])
  return(res)
}
shannon_list <- c("agar", "list")
factor_list <- readLines("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/PERMANOVA/subsample/factor_list.txt")
param_grid <- expand_grid(
  shannon = shannon_list,
  factor0 = factor_list
)
result <- pmap_dfr(param_grid, get_permanova_res)
write.table(result, "sensitive_analysis/res_permanova_macrofungi_check.tsv", sep = "\t", row.names = F, quote = F)


# --- Alpha-diversity effect-size sensitivity analysis for macrofungal filtering ---
div <- read.csv("macrofungi_check/macrofungi_filter_div_sample.csv", row.names = 1)
meta <- read.csv("/storage/zhengjushengLab/huangkailang/temp/fungi_ITS/data_for_analysis/ID_meta_final.csv", header = T, row.names = 1)

div_meta <- div %>% left_join(meta, by = "ID")
get_eta_function <- function(shannon_in, var){
  formula_lm <- as.formula(paste0(shannon_in, " ~ ", var))
  lm_model <- lm(formula_lm, data = div_meta, na.action = na.omit)

  gl <- broom::glance(lm_model)
  eta <- effectsize::eta_squared(lm_model, partial = FALSE)

  df_kruskal <- div_meta %>%
    dplyr::select(!!sym(shannon_in), !!sym(var)) %>%
    na.omit()

  result <- data.frame(
    diversity_index = shannon_in,
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
shannon_list <- c("Shannon_micro_dual", "Shannon_micro_list", "Shannon_micro_agar")
combinations <- tidyr::expand_grid(
  shannon = shannon_list,
  var = vars_to_test
)
result <- purrr::map2_dfr(
  combinations$shannon,
  combinations$var,
  get_eta_function
)
write.csv(result, "sensitive_analysis/effect_size_diversity_lm_macro.csv", row.names = F)
