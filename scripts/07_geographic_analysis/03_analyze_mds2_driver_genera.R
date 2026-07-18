# Purpose: Identify prevalent genera whose CLR abundances are most strongly associated with the major MDS axes.
# Panels: Figure S7E.

library(tidyverse)

# --- Identify genera associated with the MDS axes by Spearman correlation ---
meta <- read.csv("ID_meta_final.csv", header = T, row.names = 1)
res_mds <- read.csv("sensitive_analysis/res_mds_for_public.csv")
explained <- read.csv("sensitive_analysis/res_mds_explained_for_public.csv")
genus_clr <- read.csv("otu_gen_P0.01_filterdp_clr.csv", sep = ",", header = T, row.names = 1)
genus_clr <- t(genus_clr) %>% data.frame() %>% rownames_to_column(var = "ID")
res <- res_mds %>% left_join(genus_clr, by = "ID")

genus_list <- str_subset(names(res), "^g__")
mds_list <- c("mds1", "mds2")
get_cor_function <- function(genus0, mds0){
  df <- res %>%
    select(all_of(c("ID", mds0, genus0))) %>%
    rename_with(~c("ID", "mds", "genus"))
  cor_test_result <- cor.test(df$mds, df$genus, method = "spearman")
  res_cor <- data.frame(r_value = cor_test_result$estimate,
                        p_value = cor_test_result$p.value,
                        genus = genus0,
                        mds = mds0)
  return(res_cor)
}
result <- crossing(genus = genus_list, mds = mds_list) %>%
  pmap_dfr(~get_cor_function(genus0 = ..1, mds0 = ..2)) %>%
  arrange(mds, -r_value)

res_cor_mds2_taxa <- result %>%
  group_by(mds) %>%
  arrange(-abs(r_value)) %>%
  slice_head(n = 10) %>%
  ungroup %>%
  mutate(abs_r_value = abs(r_value)) %>%
  arrange(mds, abs(r_value)) %>%
  mutate(Genus_Unique = paste(mds, genus, sep = "__")) %>%
  mutate(Genus_Unique = factor(Genus_Unique, levels = unique(Genus_Unique))) %>%
  mutate(cor_direction = ifelse(r_value < 0, "neg", "pos")) %>%
  filter(mds == "mds2")
write.csv(res_cor_mds2_taxa, "sensitive_analysis/mds2_driver_genus.csv", row.names = F)
