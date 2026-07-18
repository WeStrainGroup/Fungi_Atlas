# Purpose: Estimate factor effect sizes for alpha diversity, community structure, genus abundance, and technical contrasts.
# Panels: Figure 3A-E. Several model-fitting sections are intended for an HPC environment.

library(phyloseq)
library(microbiome)
library(data.table)
library(Biostrings)
library(tidyverse)
library(vegan)
library(effectsize)
library(bigmds)
library(broom)
library(future)
library(future.apply)
library(clusterSim)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(iNEXT)
library(parallel)
library(rstatix)

# Read a large OTU table and restore its first column as row names.
read_otu <- function(filepath) {
  otu <- fread(filepath, header = TRUE)
  otu <- as.data.frame(otu)
  rownames(otu) <- otu[[1]]
  otu <- otu[, -1]
  return(otu)
}

# --- Assess effects of metadata factors on diversity indices ---
# Merge diversity metrics with metadata
div <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/diversity/diversity_microfungi_filterdp.csv", row.names = 1)
div <- div %>% rownames_to_column("ID") %>% dplyr::select(-PRJ)
head(div)
nrow(div)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
head(meta)
colnames(meta)
nrow(meta)

div_meta <- merge(div, meta, by = "ID", all.x = TRUE)
nrow(div_meta) #37417
colnames(div_meta)


div_list <- c("Observed", "Shannon", "Simpson")
vars_to_test <- c("DNA_extract_kit", "DNA_extract_lyticase", "DNA_extract_bead_beating",
                  "PCR_marker", "PCR_primer_pair_name", "PCR_enzyme",
                  "Seq_platform", "Seq_cycle", "Depth_microfungi_group",
                  "Sample_type", "Sex", "Age_group",
                  "Continent", "Country", "Region",
                  "PRJ")


# Run in parallel
# Important settings
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")
Sys.setenv(BLAS_NUM_THREADS = "1")

# Limit CPU usage
plan(multisession, workers = 16)

# Enumerate combinations
combo <- expand.grid(div_index = div_list, factor = vars_to_test, stringsAsFactors = FALSE)
combo

# Use a global object to avoid copying
div_meta_global <- div_meta

# -- Simple linear model (LM) --
# Parallel apply
results_list <- list()
results_list <- future_lapply(1:nrow(combo), function(i) {
  div <- combo$div_index[i]
  var <- combo$factor[i]
  cat("\n", i, " - ", div, " - ", var)

  div_meta_global[[var]] <- as.factor(div_meta_global[[var]])

  # Build model
  formula_lm <- as.formula(paste0("`", div, "` ~ ", var))
  lm_model <- lm(formula_lm, data = div_meta_global, na.action = na.omit)

  # Extract metrics
  gl <- broom::glance(lm_model)
  eta <- effectsize::eta_squared(lm_model, partial = FALSE)

  df_kruskal <- div_meta_global %>%
    dplyr::select(all_of(div), !!sym(var)) %>%
    na.omit()

  data.frame(
    div_index = div,
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
})

# Combine results
results <- do.call(rbind, results_list)
head(results)

# fdr adj
results <- results %>%
  group_by(div_index) %>%
  mutate(p_lm_adj = p.adjust(p_lm, method = "fdr"),
         p_anova_adj = p.adjust(p_anova, method = "fdr"),
         p_kruskal_adj = p.adjust(p_kruskal, method = "fdr")) %>%
  ungroup()


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/effect_size"
write.csv(results, file.path(path, paste0("effect_size_diversity_lm.csv")),  row.names = FALSE)

# --- Explain variation in community structure by factors ---
otu_gen_dist <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/data/distance/otu_gen_0.001_robait_filterdp_dist.rds")
meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)

vars_to_test <- c("DNA_extract_kit", "DNA_extract_lyticase", "DNA_extract_bead_beating",
                  "PCR_marker", "PCR_primer_pair_name", "PCR_enzyme",
                  "Seq_platform", "Seq_cycle", "Depth_microfungi_group",
                  "Sample_type", "Sex", "Age_group",
                  "Continent", "Country", "Region",
                  "PRJ")

# Align sample order to match PERMANOVA requirements
meta <- meta[match(rownames(as.matrix(otu_gen_dist)), meta$ID), ]
head(meta$ID)
head(rownames(as.matrix(otu_gen_dist)))

# Run one by one on HPC (computationally intensive)
adonis2 <- adonis2(otu_gen_dist ~ DNA_extract_kit,
                   data = meta, permutations = 999, na.action = na.omit, parallel = 2)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/effect_size"
saveRDS(adonis2, file.path(path, "adonis2_DNA_extract_kit.rds"))





# --- Assess effects of factors on single-genus abundance ---
# Load data
otu <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/genus/otu_gen_P0.01_filterdp_clr.csv", row.names = 1)
otu <- as.data.frame(t(otu))
otu[1:3, 1:3]
dim(otu)

# Remove 'Other' and 'Unknown'
otu <- otu %>% dplyr::select(-Other, -Unknown)
dim(otu)

gen_list <- colnames(otu)
gen_list #198

# Merge OTU with metadata
otu <- otu %>% rownames_to_column("ID")
head(otu)

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
head(meta)

otu_meta <- merge(otu, meta, by = "ID", all.x = TRUE)
nrow(otu_meta) #37417
colnames(otu_meta)


vars_to_test <- c("DNA_extract_kit", "DNA_extract_lyticase", "DNA_extract_bead_beating",
                  "PCR_marker", "PCR_primer_pair_name", "PCR_enzyme",
                  "Seq_platform", "Seq_cycle", "Depth_microfungi_group",
                  "Sample_type", "Sex", "Age_group",
                  "Continent", "Country", "Region",
                  "PRJ")


# Run in parallel
# Important settings
Sys.setenv(OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           BLAS_NUM_THREADS = "1")


# Limit CPU usage
plan(multisession, workers = 32)

# Enumerate combinations
combo <- expand.grid(genus = gen_list, factor = vars_to_test, stringsAsFactors = FALSE)
combo

# Use a global object to avoid copying
otu_meta_global <- otu_meta

results_list <- list()
results_list <- future_lapply(1:nrow(combo), function(i) {
  g <- combo$genus[i]
  var <- combo$factor[i]
  cat(paste0("\n[", i, "] ", g, " ~ ", var))

  otu_meta_global[[var]] <- as.factor(otu_meta_global[[var]])

  # Build model
  formula_lm <- as.formula(paste0("`", g, "` ~ ", var))
  lm_model <- lm(formula_lm, data = otu_meta_global, na.action = na.omit)
  gl <- broom::glance(lm_model)

  # Eta-squared
  eta <- effectsize::eta_squared(lm_model, partial = FALSE)

  # Kruskal–Wallis test
  df_kruskal <- otu_meta_global %>%
    dplyr::select(all_of(g), !!sym(var)) %>%
    na.omit()
  p_kruskal <- kruskal.test(formula_lm, data = df_kruskal)$p.value

  # Return results
  data.frame(
    genus = g,
    factor = var,
    n_sample = gl$df.residual + 1,
    r_squared = gl$r.squared,
    r_squared_adj = gl$adj.r.squared,
    eta_squared = eta$Eta2[eta$Parameter == var],
    p_lm = gl$p.value,
    p_anova = anova(lm_model)$`Pr(>F)`[1],
    p_kruskal = p_kruskal,
    stringsAsFactors = FALSE
  )
})

# Combine results
results <- do.call(rbind, results_list)

# fdr adj
results <- results %>%
  group_by(genus) %>%
  mutate(p_lm_adj = p.adjust(p_lm, method = "fdr"),
         p_anova_adj = p.adjust(p_anova, method = "fdr"),
         p_kruskal_adj = p.adjust(p_kruskal, method = "fdr")) %>%
  ungroup()


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/effect_size"
write.csv(results, file.path(path, paste0("effect_size_genus_lm.csv")), row.names = FALSE)




# --- Differential analysis: effects of various factors on OTUs ---
otu <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/genus/otu_gen_P0.01_filterdp_clr.csv")
otu <- t(otu) %>%
    as.data.frame() %>%
    rownames_to_column("ID")

otu[1:3, 1:3]

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
head(meta)
colnames(meta)
nrow(meta)

meta_PRJ <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/PRJ_meta_tech.csv", row.names = 1)




# -- Differential taxa associated with PCR enzymes --
# Filter enzymes that meet coverage criteria
# Helper to summarize factor coverage
summarise_factor <- function(meta, meta_PRJ, factor) {

  # Number of projects
  df_proj <- as.data.frame(table(meta_PRJ[[factor]]))
  colnames(df_proj) <- c(factor, "n_PRJ")

  # Number of samples
  df_sample <- as.data.frame(table(meta[[factor]]))
  colnames(df_sample) <- c(factor, "n_sample")

  # Number of continents (based on 'Continent')
  df_continent <- meta %>%
    filter(!is.na(Continent), !is.na(.data[[factor]])) %>%
    group_by(.data[[factor]], Continent) %>%
    summarise(n_sample = n(), .groups = "drop") %>%
    group_by(.data[[factor]]) %>%
    summarise(n_continent = n_distinct(Continent), .groups = "drop")

  # Merge summaries
  df_merge1 <- merge(df_proj, df_sample, by = factor, all = F)
  df_merge <- merge(df_merge1, df_continent, by = factor, all = F)

  return(df_merge)
}

# Summaries
df_enzyme <- summarise_factor(meta, meta_PRJ, "PCR_enzyme")
head(df_enzyme)

df_enzyme_filter <- df_enzyme %>%
    filter(n_PRJ >= 5 & n_sample >= 500 & n_continent >= 3)
df_enzyme_filter # 2 enzymes

# Subset samples with selected enzymes
otu_meta_filter <- otu_meta %>%
    filter(PCR_enzyme %in% df_enzyme_filter$PCR_enzyme)

unique(otu_meta_filter$PCR_enzyme)
nrow(otu_meta_filter) # 2364

# Cast to factor
otu_meta_filter[["PCR_enzyme"]] <- as.factor(otu_meta_filter[["PCR_enzyme"]])

# Linear mixed model (LMM)
results_list <- list()
results_list <- lapply(taxa_list, function(taxon) {
  # Build model
  formula <- as.formula(paste(taxon, "~ PCR_enzyme + (1|PRJ)"))
  fit <- lmer(formula, data = otu_meta_filter)

  # Tidy results
  tidy_fit <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
  tidy_fit <- tidy_fit %>% filter(term == "PCR_enzymePhusion High-Fidelity")
  tidy_fit$taxon <- taxon   # Add taxon name

  tidy_fit # Return
})

res_lmm <- bind_rows(results_list)
res_lmm <- res_lmm %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  arrange(p_adj)
head(res_lmm)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/differential_analysis"
write.csv(res_lmm, file.path(path, "lmm_PCR_enzyme.csv"))



# -- Differential taxa associated with PCR primer pairs --
# Summaries
df_primer <- summarise_factor(meta, meta_PRJ, "PCR_primer_pair_name")
head(df_enzyme)


df_primer_filter <- df_primer %>%
    filter(n_PRJ >= 10 & n_sample >= 1000 & n_continent >= 3)
df_primer_filter # 2 primer pairs

otu_meta_filter <- otu_meta %>%
    filter(PCR_primer_pair_name %in% df_primer_filter$PCR_primer)

nrow(otu_meta)
nrow(otu_meta_filter) # 33581

# Cast to factor
otu_meta_filter[["PCR_primer_pair_name"]] <- as.factor(otu_meta_filter[["PCR_primer_pair_name"]])

results_list <- list()
results_list <- lapply(taxa_list, function(taxon) {
  # Build model
  formula <- as.formula(paste(taxon, "~ PCR_primer_pair_name + (1|PRJ)"))
  fit <- lmer(formula, data = otu_meta_filter)

  # Tidy results
  tidy_fit <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
  tidy_fit <- tidy_fit %>% filter(term == "PCR_primer_pair_nameITS3_ITS4")
  tidy_fit$taxon <- taxon   # Add taxon name

  tidy_fit # Return
})


res_lmm <- bind_rows(results_list)
res_lmm <- res_lmm %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  arrange(p_adj)
head(res_lmm)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/differential_analysis"
write.csv(res_lmm, file.path(path, "lmm_PCR_primer.csv"))



# -- Sequencing platform differences --
df_platform <- summarise_factor(meta, meta_PRJ, "Seq_platform")
head(df_platform)

# Directly select platforms
otu_meta_filter <- otu_meta %>%
    filter(Seq_platform %in% c("MiSeq", "NovaSeq"))

nrow(otu_meta)
nrow(otu_meta_filter) #36542

# Cast to factor
otu_meta_filter[["Seq_platform"]] <- as.factor(otu_meta_filter[["Seq_platform"]])


results_list <- list()
results_list <- lapply(taxa_list, function(taxon) {
  # Build model
  formula <- as.formula(paste(taxon, "~ Seq_platform + (1|PRJ)"))
  fit <- lmer(formula, data = otu_meta_filter)

  # Tidy results
  tidy_fit <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE)
  tidy_fit <- tidy_fit %>% filter(term == "Seq_platformNovaSeq")   # Keep the 'Seq_platformNovaSeq' term
  tidy_fit$taxon <- taxon   # Add taxon name

  tidy_fit # Return
})


res_lmm <- bind_rows(results_list)
res_lmm <- res_lmm %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  arrange(p_adj)
head(res_lmm)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/differential_analysis"
write.csv(res_lmm, file.path(path, "lmm_Seq_platform.csv"))
