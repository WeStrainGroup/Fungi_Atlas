# Purpose: Test cohort-level diversity, run ANCOM-BC2, and evaluate disease classifiers at genus and SH levels.
# Panels: Figure 6A-C; Figure S9A-C.

library(tidyverse)
library(patchwork)
library(vegan)
library(ggplot2)
library(hugmycoa)
library(phyloseq)

setwd("/data/huangkailang/project/fungi_ITS/case_control_metadata_clean")
sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)



### alpha diversity wilcox test
ps <- load_hugmycoa(level = "Genus", filterdepth = TRUE)
ps_otu <- otu_table(ps)
ps_tax <- tax_table(ps)
ps_meta <- sample_data(ps) %>%
    data.frame() %>%
    mutate(project = str_split_i(PRJ, "_", 1))  %>%
    relocate(project, .after = PRJ)

alpha_diversity <- ps_meta %>%
    select(project, ID, Shannon, Simpson, Observed)

metadata <- read.csv("table/metadata_with_filterdepth.csv")
cohort_list <- cohort_order
alpha_diversity_function <- function(cohort0) {
    id0 <- metadata %>% filter(cohort == cohort0)  %>% pull(sample)
    meta0 <- metadata %>% filter(cohort == cohort0)
    alpha0 <- alpha_diversity %>%
        filter(ID %in% id0) %>%
        left_join(meta0, by = c("ID" = "sample")) %>%
        mutate(group = factor(group, levels = c("Case", "Control")))

    wilcox_result1 <- wilcox.test(Shannon ~ group, data = alpha0, conf.int = TRUE)
    wilcox_result2 <- wilcox.test(Simpson ~ group, data = alpha0, conf.int = TRUE)
    wilcox_result3 <- wilcox.test(Observed ~ group, data = alpha0, conf.int = TRUE)

    forest_df <- data.frame(
        alpha_index = c("Shannon", "Simpson", "Observed"),
        Estimate = c(wilcox_result1$estimate, wilcox_result2$estimate, wilcox_result3$estimate),
        CI_low = c(wilcox_result1$conf.int[1], wilcox_result2$conf.int[1], wilcox_result3$conf.int[1]),
        CI_high = c(wilcox_result1$conf.int[2], wilcox_result2$conf.int[2], wilcox_result3$conf.int[2]),
        P_value = c(wilcox_result1$p.value, wilcox_result2$p.value, wilcox_result3$p.value)
        ) %>%
        mutate(cohort = cohort0)
    return(forest_df)
}
result_alpha_all <- map_dfr(cohort_list, alpha_diversity_function)
result_alpha_all <- result_alpha_all %>%
    mutate(significance = case_when(
      P_value < 0.001 ~ "***",
      P_value < 0.01 ~ "**",
      P_value < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    mutate(alpha_index = factor(alpha_index, levels = c("Shannon", "Simpson", "Observed"))) %>%
    arrange(alpha_index, Estimate)
head(result_alpha_all)
write.csv(result_alpha_all, "table/res_alpha_diversity_wilcox_test.csv", row.names = FALSE)



### beta diversity adonis test (genus level)
ps <- load_hugmycoa(level = "Genus", filterdepth = TRUE)
ps_otu <- otu_table(ps)
metadata <- read.csv("table/metadata_with_filterdepth.csv")
cohort_list <- cohort_order

beta_diversity_function <- function(cohort0, distance_method) {
  id0 <- metadata %>% filter(cohort == cohort0)  %>% pull(sample) %>% intersect(sample_names(ps_otu))
  meta0 <- metadata %>% filter(cohort == cohort0) %>% filter(sample %in% id0)

  abun0 <- ps_otu[, id0, drop = FALSE] %>%
    data.frame() %>%
    filter(rowSums(. > 0) > 0 )
  detect_taxa <- rownames(abun0)

  dist_robust <- vegdist(t(abun0), method = distance_method)
  pcoa_res <- cmdscale(dist_robust, k = 10, eig = TRUE)

  # explained
  pcoa_eig <- pcoa_res$eig
  pcoa_eig <- round((pcoa_eig / sum(pcoa_eig)) * 100, 1)

  pcoa_point <- as.data.frame(pcoa_res$points) %>%
    rename_with(~ paste0("PCoA", 1:10)) %>%
    rownames_to_column("sample") %>%
    select(sample, PCoA1, PCoA2) %>%
    left_join(meta0, by = "sample") %>%
    mutate(pcoa_eig1 = pcoa_eig[1], pcoa_eig2 = pcoa_eig[2] ) %>%
    relocate(cohort, .before = sample)

  # PERMANOVA
  pcoa_point <- pcoa_point[match(rownames(as.matrix(dist_robust)), pcoa_point$sample), ]
  head(pcoa_point$sample);length(pcoa_point$sample);
  head(rownames(as.matrix(dist_robust))); length(rownames(as.matrix(dist_robust)))
  set.seed(123)
  adonis2_res <- adonis2(dist_robust ~ group,
                         data = pcoa_point,
                         permutations = 999)
  results <- pcoa_point %>%
    mutate(adonis_R2 = adonis2_res$R2[1], adonis_p_value = adonis2_res$`Pr(>F)`[1])
  return(results)
}

result_beta_all_robust_genus <- map_dfr(cohort_list, beta_diversity_function, distance_method = "robust.aitchison")
write.csv(result_beta_all_robust_genus, "table/res_beta_diversity_pcoa_genus.csv", row.names = FALSE)

result_adonis_robust_genus <- result_beta_all_robust_genus %>%
    select(cohort, adonis_R2, adonis_p_value) %>%
    unique()  %>%
    mutate(significance = case_when(
      adonis_p_value < 0.001 ~ "***",
      adonis_p_value < 0.01 ~ "**",
      adonis_p_value < 0.05 ~ "*",
      TRUE ~ ""
    ))  %>%
    mutate(cohort = factor(cohort, levels = cohort_order))  %>%
    arrange(cohort)
write.csv(result_adonis_robust_genus, "table/res_beta_diversity_adonis_test_robust_genus.csv", row.names = FALSE)



### beta diversity adonis test (SH level)
ps_clu <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Main/genetic_distance/clustered/ps_clu_fdp.rds")
ps_otu <- otu_table(ps_clu) %>% t()

metadata <- read.csv("table/metadata_with_filterdepth.csv")
sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)
cohort_list <- cohort_order

beta_diversity_function <- function(cohort0, distance_method) {
  id0 <- metadata %>% filter(cohort == cohort0)  %>% pull(sample) %>% intersect(sample_names(ps_otu))
  meta0 <- metadata %>% filter(cohort == cohort0) %>% filter(sample %in% id0)

  abun0 <- ps_otu[, id0, drop = FALSE] %>%
    data.frame() %>%
    filter(rowSums(. > 0) > 0 )
  detect_taxa <- rownames(abun0)

  dist_robust <- vegdist(t(abun0), method = distance_method) # row: samples，column: feature
  pcoa_res <- cmdscale(dist_robust, k = 10, eig = TRUE)

  # explained
  pcoa_eig <- pcoa_res$eig
  pcoa_eig <- round((pcoa_eig / sum(pcoa_eig)) * 100, 1)

  pcoa_point <- as.data.frame(pcoa_res$points) %>%
    rename_with(~ paste0("PCoA", 1:10)) %>%
    rownames_to_column("sample") %>%
    select(sample, PCoA1, PCoA2) %>%
    left_join(meta0, by = "sample") %>%
    mutate(pcoa_eig1 = pcoa_eig[1], pcoa_eig2 = pcoa_eig[2] ) %>%
    relocate(cohort, .before = sample)

  # PERMANOVA
  pcoa_point <- pcoa_point[match(rownames(as.matrix(dist_robust)), pcoa_point$sample), ]
  head(pcoa_point$sample);length(pcoa_point$sample);
  head(rownames(as.matrix(dist_robust))); length(rownames(as.matrix(dist_robust)))
  set.seed(123)
  adonis2_res <- adonis2(dist_robust ~ group,
                         data = pcoa_point,
                         permutations = 999)
  results <- pcoa_point %>%
    mutate(adonis_R2 = adonis2_res$R2[1], adonis_p_value = adonis2_res$`Pr(>F)`[1])
  return(results)
}
result_beta_all_robust_sh <- map_dfr(cohort_list, beta_diversity_function, distance_method = "robust.aitchison")
write.csv(result_beta_all_robust_sh, "table/res_beta_diversity_pcoa_sh.csv", row.names = FALSE)

result_adonis_robust_sh <- result_beta_all_robust_sh %>%
    select(cohort, adonis_R2, adonis_p_value) %>%
    unique()  %>%
    mutate(significance = case_when(
      adonis_p_value < 0.001 ~ "***",
      adonis_p_value < 0.01 ~ "**",
      adonis_p_value < 0.05 ~ "*",
      TRUE ~ ""
    ))  %>%
    mutate(cohort = factor(cohort, levels = cohort_order))  %>%
    arrange(cohort)
write.csv(result_adonis_robust_sh, "table/res_beta_diversity_adonis_test_robust_sh.csv", row.names = FALSE)



### relative abundance differential analysis （genus level）
library(ANCOMBC)
ps <- load_hugmycoa(level = "Genus", filterdepth = TRUE)
ps_otu <- otu_table(ps)
ps_otu[1:5, 1:5]

metadata <- read.csv("table/metadata_with_filterdepth.csv")
sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)
cohort_list <- cohort_order

ancombc2_function <- function(cohort0) {
  id0 <- metadata %>% filter(cohort == cohort0)  %>% pull(sample) %>% intersect(sample_names(ps_otu))
  meta0 <- metadata %>% filter(cohort == cohort0) %>%
    filter(sample %in% id0) %>%
    column_to_rownames("sample")

  abun0 <- ps_otu[, id0, drop = FALSE] %>%
    data.frame() %>%
    filter(rowSums(. > 0) > 0 * ncol(.))
  detect_taxa <- rownames(abun0)

  meta0$group <- factor(meta0$group, levels = c("Control", "Case"))
  res0 <- ancombc2(
    data = abun0,
    meta_data = meta0,
    fix_formula = "group",
    group = "group",
    p_adj_method = "BH",
    prv_cut = 0.1, # prevalence filter
    lib_cut = 0, # sequence depth filter
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05 # FDR threshold
  )

  result_anc0 <- res0$res %>%
    dplyr::select(taxon, lfc_groupCase, p_groupCase, q_groupCase, diff_groupCase, passed_ss_groupCase) %>%
    rename(feature_id = taxon, log2FC = lfc_groupCase, p_value = p_groupCase, q_value = q_groupCase, diff = diff_groupCase, passed_ss = passed_ss_groupCase) %>%
    bind_rows(data.frame(feature_id = setdiff(detect_taxa, res0$res$taxon))) %>%
    mutate(cohort = cohort0) %>%
    relocate(cohort, .before = feature_id)
  return(result_anc0)
}
res_diff_genus <- map_dfr(cohort_list, ancombc2_function) %>%
  filter(!feature_id %in% c("Unknown"))
write.csv(res_diff_genus, "table/ancombc2_results_genus.csv", row.names = F)



### relative abundance differential analysis （SH level）
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(ANCOMBC)

ps_clu <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Main/genetic_distance/clustered/ps_clu_fdp.rds")
ps_otu <- otu_table(ps_clu) %>% t()
ps_otu[1:5, 1:5]

ps_tax <- tax_table(ps_clu)
sh_label_table <- ps_tax %>%
  data.frame() %>%
  rownames_to_column(var = "SH") %>%
  mutate(sh_label = case_when(!is.na(Species) ~ paste(str_remove_all(Genus, "g__"), str_remove_all(Species, "s__"), SH, sep = " "),
                              !is.na(Genus) ~ paste(Genus, "sp.", SH, sep = " "),
                              !is.na(Family) ~ paste(Family, "sp.", SH, sep = " "),
                              !is.na(Order) ~ paste(Order, "sp.", SH, sep = " "),
                              !is.na(Class) ~ paste(Class, "sp.", SH, sep = " "),
                              !is.na(Phylum) ~ paste(Phylum, "sp.", SH, sep = " "),
                              T ~ paste0("unclassified ", SH)) %>%
           str_remove_all("p__|c__|o__|f__|g__|s__")
        )
write.csv(sh_label_table, "table/sh_label_table.csv", row.names = F)

metadata <- read.csv("table/metadata_with_filterdepth.csv")
sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)
cohort_list <- cohort_order

ancombc2_function <- function(cohort0) {
  id0 <- metadata %>% filter(cohort == cohort0)  %>% pull(sample) %>% intersect(sample_names(ps_otu))
  meta0 <- metadata %>% filter(cohort == cohort0) %>%
    filter(sample %in% id0) %>%
    column_to_rownames("sample")

  abun0 <- ps_otu[, id0, drop = FALSE] %>%
    data.frame() %>%
    filter(rowSums(. > 0) > 0 * ncol(.))
  detect_taxa <- rownames(abun0)

  meta0$group <- factor(meta0$group, levels = c("Control", "Case"))
  res0 <- ancombc2(
    data = abun0,
    meta_data = meta0,
    fix_formula = "group",
    group = "group",
    p_adj_method = "BH",
    prv_cut = 0.1,
    lib_cut = 0,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05
  )
  head(res0)

  result_anc0 <- res0$res %>%
    dplyr::select(taxon, lfc_groupCase, p_groupCase, q_groupCase, diff_groupCase, passed_ss_groupCase) %>%
    rename(feature_id = taxon, log2FC = lfc_groupCase, p_value = p_groupCase, q_value = q_groupCase, diff = diff_groupCase, passed_ss = passed_ss_groupCase) %>%
    bind_rows(data.frame(feature_id = setdiff(detect_taxa, res0$res$taxon))) %>%
    mutate(cohort = cohort0) %>%
    relocate(cohort, .before = feature_id)
  head(result_anc0)
  return(result_anc0)
}
res_diff_sh <- map_dfr(cohort_list, ancombc2_function) %>%
  filter(!feature_id %in% c("Unknown"))
write.csv(res_diff_sh, "table/ancombc2_results_SH.csv", row.names = F)



### prediction based on different genus
library(tidyverse)
library(randomForest)
library(glmnet)
library(pROC)
library(hugmycoa)
library(phyloseq)

ps <- load_hugmycoa(level = "Genus", filterdepth = TRUE)
ps_otu <- otu_table(ps)
metadata <- read.csv("table/metadata_with_filterdepth.csv")
sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)
cohort_list <- cohort_order

# different taxa
res_diff_genus <- read.csv("table/ancombc2_results_genus.csv") %>%
  filter(abs(log2FC) > 1 & q_value < 0.05)

# random forest, svm, lasso
library(e1071)
library(purrr)
library(dplyr)

get_clean_data <- function(cohort0) {
  # 1. get diff taxa and samples
  id0 <- metadata %>% filter(cohort == cohort0)  %>% pull(sample) %>% intersect(sample_names(ps_otu))
  meta0 <- metadata %>% filter(cohort == cohort0) %>% filter(sample %in% id0)
  feature0 <- res_diff_genus %>% filter(cohort == cohort0) %>% pull(feature_id)

  if (length(feature0) == 0) {return(NULL)}

  otu_all_raw <- as.data.frame(t(ps_otu[, id0] %>% data.frame() %>% filter(rowSums(. > 0) > 0.05 * ncol(.))))
  otu_all_prop <- t(apply(otu_all_raw, 1, function(x) (x + 1) / sum(x + 1)))
  otu_all_clr <- t(apply(otu_all_prop, 1, function(x) { log_x <- log(x); return(log_x - mean(log_x)) })) %>% as.data.frame()

  meta1 <- meta0 %>% select(sample, group) %>% column_to_rownames(var = "sample") %>% mutate(group = factor(group, levels = c("Control", "Case")))
  y <- factor(meta1$group, levels = c("Control", "Case"))
  names(y) <- rownames(meta1)

  feature0_clean <- intersect(feature0, colnames(otu_all_clr))
  X <- otu_all_clr %>% select(all_of(feature0_clean)) %>% filter(rownames(.) %in% id0)
  X <- X[names(y), , drop = FALSE]

  return(list(X = X, y = y))
}


run_loocv <- function(X, y, method) {
  n <- length(y)
  p <- ncol(X)

  probs <- sapply(1:n, function(i) {
    X_tr <- scale(X[-i, , drop=FALSE])
    X_te <- scale(X[i, , drop=FALSE], center=attr(X_tr,"scaled:center"), scale=attr(X_tr,"scaled:scale"))
    X_te[is.na(X_te)] <- 0
    set.seed(i)

    # 1. random forest
    if(method == "rf") {
      return(predict(randomForest(x = X_tr, y = y[-i], ntree=300), newdata = X_te, type = "prob")[,"Case"])
    }

    # 2. Lasso / GLM
    if(method == "lasso") {
      if (p <= 1) {
        train_df <- data.frame(Feature = as.numeric(X_tr), group = y[-i], row.names = NULL)
        test_df  <- data.frame(Feature = as.numeric(X_te), row.names = NULL)
        return(predict(glm(group ~ Feature, data = train_df, family = binomial), newdata = test_df, type = "response"))
      } else {
        return(as.numeric(predict(cv.glmnet(X_tr, y[-i], family="binomial", alpha=1), newx = X_te, s="lambda.min", type="response")))
      }
    }
  })

  # results
  rf_roc <- roc(y, probs, levels = c("Control", "Case"), direction = "<", quiet = TRUE)
  return(data.frame(auc = round(auc(rf_roc), 3), low = round(ci.auc(rf_roc)[1], 3), high = round(ci.auc(rf_roc)[3], 3)))
}


set.seed(123)
final_comparison <- map_dfr(cohort_list, function(coh) {
  dat <- get_clean_data(coh)
  if(is.null(dat)) return(data.frame(cohort = rep(coh, 2), model = c("rf", "lasso"), n_features = rep(0, 2), auc = rep(0, 2), low = rep(0, 2), high = rep(0, 2)))  #svm
  map_dfr(c("rf", "lasso"), function(m) {
    res <- run_loocv(dat$X, dat$y, m)
    data.frame(cohort = coh, model = m, n_features = ncol(dat$X), res)
  })
}) %>% arrange(cohort, -auc)
res_prediction_genus <- final_comparison
head(res_prediction_genus)
write.csv(res_prediction_genus, "table/res_prediction_genus.csv", row.names = F)



### prediction based on different SH
library(hugmycoa)
library(phyloseq)

ps <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Main/genetic_distance/clustered/ps_clu_fdp.rds")
ps_otu <- otu_table(ps) %>% t()

metadata <- read.csv("table/metadata_with_filterdepth.csv")
sample_number <- read.csv("table/sample_number_filterdepth.csv")
cohort_order <- unique(sample_number$cohort)
cohort_list <- cohort_order

# different taxa
res_diff_sh <- read.csv("table/ancombc2_results_SH.csv") %>%
  filter(abs(log2FC) > 1 & q_value < 0.05)

# random forest, svm, lasso
library(e1071)
library(purrr)
library(dplyr)

get_clean_data <- function(cohort0) {
  # get diff taxa and samples
  id0 <- metadata %>% filter(cohort == cohort0)  %>% pull(sample) %>% intersect(sample_names(ps_otu))
  meta0 <- metadata %>% filter(cohort == cohort0) %>% filter(sample %in% id0)
  feature0 <- res_diff_sh %>% filter(cohort == cohort0) %>% pull(feature_id)

  if (length(feature0) == 0) {return(NULL)}

  otu_all_raw <- as.data.frame(t(ps_otu[, id0] %>% data.frame() %>% filter(rowSums(. > 0) > 0.05 * ncol(.))))
  otu_all_prop <- t(apply(otu_all_raw, 1, function(x) (x + 1) / sum(x + 1)))
  otu_all_clr <- t(apply(otu_all_prop, 1, function(x) { log_x <- log(x); return(log_x - mean(log_x)) })) %>% as.data.frame()

  meta1 <- meta0 %>% select(sample, group) %>% column_to_rownames(var = "sample") %>% mutate(group = factor(group, levels = c("Control", "Case")))
  y <- factor(meta1$group, levels = c("Control", "Case"))
  names(y) <- rownames(meta1)

  feature0_clean <- intersect(feature0, colnames(otu_all_clr))
  X <- otu_all_clr %>% select(all_of(feature0_clean)) %>% filter(rownames(.) %in% id0)
  X <- X[names(y), , drop = FALSE]

  return(list(X = X, y = y))
}

run_loocv <- function(X, y, method) {
  n <- length(y)
  p <- ncol(X)

  probs <- sapply(1:n, function(i) {
    X_tr <- scale(X[-i, , drop=FALSE])
    X_te <- scale(X[i, , drop=FALSE], center=attr(X_tr,"scaled:center"), scale=attr(X_tr,"scaled:scale"))
    X_te[is.na(X_te)] <- 0
    set.seed(i)

    # 1. random forest
    if(method == "rf") {
      return(predict(randomForest(x = X_tr, y = y[-i], ntree=300), newdata = X_te, type = "prob")[,"Case"])
    }

    # 2. Lasso / GLM
    if(method == "lasso") {
      if (p <= 1) {
        train_df <- data.frame(Feature = as.numeric(X_tr), group = y[-i], row.names = NULL)
        test_df  <- data.frame(Feature = as.numeric(X_te), row.names = NULL)
        return(predict(glm(group ~ Feature, data = train_df, family = binomial), newdata = test_df, type = "response"))
      } else {
        return(as.numeric(predict(cv.glmnet(X_tr, y[-i], family="binomial", alpha=1), newx = X_te, s="lambda.min", type="response")))
      }
    }
  })

  # results
  rf_roc <- roc(y, probs, levels = c("Control", "Case"), direction = "<", quiet = TRUE)
  return(data.frame(auc = round(auc(rf_roc), 3), low = round(ci.auc(rf_roc)[1], 3), high = round(ci.auc(rf_roc)[3], 3)))
}

set.seed(123)
final_comparison <- map_dfr(cohort_list, function(coh) {
  dat <- get_clean_data(coh)
  if(is.null(dat)) return(data.frame(cohort = rep(coh, 2), model = c("rf", "lasso"), n_features = rep(0, 2), auc = rep(0, 2), low = rep(0, 2), high = rep(0, 2)))

  map_dfr(c("rf", "lasso"), function(m) {
    res <- run_loocv(dat$X, dat$y, m)
    data.frame(cohort = coh, model = m, n_features = ncol(dat$X), res)
  })
}) %>% arrange(cohort, -auc)
res_prediction_SH <- final_comparison
write.csv(res_prediction_SH, "table/res_prediction_SH.csv", row.names = F)
