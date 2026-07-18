# Purpose: Compare edible-mushroom-list, Agaricomycetes, and dual macrofungal filters and their diversity effects.
# Panels: Figure S4.


library(phyloseq)
library(microbiome)
library(data.table)
library(Biostrings)
library(tidyverse)
library(ggpubr)
library(vegan)




project_root1 <- "/data/wangxinyu/Fungi_Atlas/Data/Public/Annotation/SE"
project_root2 <- "/data/wangxinyu/Fungi_Atlas/Data/Public/Annotation/PE"
project_root3 <- "/data/wangxinyu/Fungi_Atlas/Data/WeGut"

project1 <- list.dirs(project_root1, full.names = TRUE, recursive = FALSE)
project2 <- list.dirs(project_root2, full.names = TRUE, recursive = FALSE)
project3 <- list.dirs(project_root3, full.names = TRUE, recursive = FALSE)

projects <- c(project1, project2, project3)
projects
length(projects) # 83


ps_merged <- NULL
not_merged_projects <- character()

for (proj in projects) {
  project_name <- basename(proj)
  cat("Processing:", project_name, "\n")

  seqtab_file <- file.path(proj, "3_phyloseq/all", paste0(project_name, "_otu_all.csv"))
  taxa_file   <- file.path(proj, "3_phyloseq/all", paste0(project_name, "_taxa_all.csv"))
  refseq_file <- file.path(proj, "3_phyloseq/all", paste0(project_name, "_refseq_all.fasta"))

  if (file.exists(seqtab_file) && file.exists(taxa_file)) {
    seqtab <- fread(seqtab_file, header = T)
    seqtab <- as.data.frame(seqtab)
    rownames(seqtab) <- seqtab[[1]]
    seqtab <- seqtab[, -1]

    taxa <- fread(taxa_file, header = T)
    taxa <- as.data.frame(taxa)
    rownames(taxa) <- taxa[[1]]
    taxa <- taxa[, -1]

    samdf <- data.frame(cohort = project_name,
                        id = rownames(seqtab),
                        row.names = rownames(seqtab))

    refseq <- readDNAStringSet(refseq_file, format = "fasta")

    ps <- phyloseq(otu_table(as.matrix(seqtab), taxa_are_rows = FALSE),
                   tax_table(as.matrix(taxa)),
                   sample_data(samdf),
                   refseq)


    if (is.null(ps_merged)) {
      ps_merged <- ps
    } else {
      ps_merged <- merge_phyloseq(ps_merged, ps)
      rm(ps)
      gc()
    }
  } else {
    cat("No file found: ", project_name, "\n")
    not_merged_projects <- c(not_merged_projects, project_name)
  }
}


not_merged_projects
ps_merged # 174175 taxa and 46598 samples

path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
saveRDS(ps_merged, file = file.path(path, "ps_all.rds"))


div_all <- estimate_richness(ps_merged, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
head(div_all)

path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
write.csv(div_all, file.path(path, "div_all.csv"))





ps_all <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check/ps_all.rds")
ps_fungi <- subset_taxa(ps_all, Kingdom == "k__Fungi")

# Check
ntaxa(ps_fungi)
get_taxa_unique(ps_fungi, "Phylum")




mushroom_list <- read.csv("/data/wangxinyu/Tools_Dev/MycoGAP/mycogap/ref/mushroom_genus.csv")
mushroom_list <- mushroom_list$Genus
mushroom_list <- paste0("g__", mushroom_list)
head(mushroom_list)

# Find ASVs whose genus name matches entries in the mushroom database
taxa_table <- as.data.frame(ps_fungi@tax_table)
taxa_table_mushroom1 <- taxa_table %>% filter(Genus %in% mushroom_list)
head(taxa_table_mushroom1)
nrow(taxa_table_mushroom1) # 13077


# Find ASVs in taxa_table where class is Agaricomycetes
taxa_table_mushroom2 <- taxa_table %>% filter(Class == "c__Agaricomycetes")
head(taxa_table_mushroom2)
nrow(taxa_table_mushroom2) # 16816

# Remove identified ASVs
ASV_to_drop <- union(rownames(taxa_table_mushroom1), rownames(taxa_table_mushroom2))
head(ASV_to_drop)
length(ASV_to_drop) # 18988



ps_macrofungi_dual <- prune_taxa((taxa_names(ps_fungi) %in% ASV_to_drop), ps_fungi)
ps_macrofungi_list <- prune_taxa((taxa_names(ps_fungi) %in% rownames(taxa_table_mushroom1)), ps_fungi)
ps_macrofungi_agar <- prune_taxa((taxa_names(ps_fungi) %in% rownames(taxa_table_mushroom2)), ps_fungi)


ps_macrofungi_dual <- filter_taxa(ps_macrofungi_dual, function(otu) sum(otu) > 1, TRUE)
ps_macrofungi_list <- filter_taxa(ps_macrofungi_list, function(otu) sum(otu) > 1, TRUE)
ps_macrofungi_agar <- filter_taxa(ps_macrofungi_agar, function(otu) sum(otu) > 1, TRUE)



path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
saveRDS(ps_macrofungi_dual, file = file.path(path, "ps_macrofungi_dual.rds"))
saveRDS(ps_macrofungi_list, file = file.path(path, "ps_macrofungi_list.rds"))
saveRDS(ps_macrofungi_agar, file = file.path(path, "ps_macrofungi_agar.rds"))

ps_macrofungi_dual <- readRDS(file.path(path, "ps_macrofungi_dual.rds"))
ps_macrofungi_list <- readRDS(file.path(path, "ps_macrofungi_list.rds"))
ps_macrofungi_agar <- readRDS(file.path(path, "ps_macrofungi_agar.rds"))



depth_macro_dual <- readcount(ps_macrofungi_dual) %>% as.data.frame() %>% rename("depth_macro_dual" = ".")
depth_macro_list <- readcount(ps_macrofungi_list) %>% as.data.frame() %>% rename("depth_macro_list" = ".")
depth_macro_agar <- readcount(ps_macrofungi_agar) %>% as.data.frame() %>% rename("depth_macro_agar" = ".")
head(depth_macro_dual)


obs_macro_dual <- estimate_richness(ps_macrofungi_dual, split = TRUE, measures = c("Observed")) %>% rename("obs_macro_dual" = "Observed")
obs_macro_list <- estimate_richness(ps_macrofungi_list, split = TRUE, measures = c("Observed")) %>% rename("obs_macro_list" = "Observed")
obs_macro_agar <- estimate_richness(ps_macrofungi_agar, split = TRUE, measures = c("Observed")) %>% rename("obs_macro_agar" = "Observed")
head(obs_macro_dual)


res_list <- list(
  depth_macro_dual,
  depth_macro_list,
  depth_macro_agar,
  obs_macro_dual,
  obs_macro_list,
  obs_macro_agar
)

depth_obs <- res_list %>%
  map(~ rownames_to_column(.x, "ID")) %>%
  reduce(full_join, by = "ID")

head(depth_obs)
dim(depth_obs)



depth_all <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Main/data/depth/depth_all.csv", row.names = 1)

div_all <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check/div_all.csv", row.names = 1)
obs_all <- div_all %>%
  select(Observed) %>%
  rename("obs_all" = "Observed") %>%
  rownames_to_column("ID")

meta <- as.matrix(sample_data(ps_merged)) %>%
  as.data.frame() %>%
  select(cohort) %>%
  rownames_to_column("ID")

head(depth_all)
nrow(depth_all)

head(obs_all)
nrow(obs_all)

head(meta)
nrow(meta)


meta_all <- merge(depth_all, meta, by = "ID")
meta_all <- merge(meta_all, obs_all, by = "ID")
meta_all <- meta_all %>% select(ID, cohort, depth_all, obs_all)
head(meta_all)


meta_all_filter <- meta_all %>% filter(depth_all >= 10000)
head(meta_all_filter)
nrow(meta_all_filter) # 42366


df <- merge(meta_all_filter, depth_obs, by = "ID", all.x = T)
head(df)
nrow(df) # 42366

path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
write.csv(df, file.path(path, "macrofungi_filter_depth_obs.csv"))



df_depth_prj <- df %>%
  mutate(
    dual = depth_macro_dual / depth_all * 100,
    list = depth_macro_list / depth_all * 100,
    agar = depth_macro_agar / depth_all * 100
  ) %>%
  group_by(cohort) %>%
  summarise(
    dual = mean(dual, na.rm = TRUE),
    list = mean(list, na.rm = TRUE),
    agar = mean(agar, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(dual, list, agar),
    names_to = "strategy",
    values_to = "percent"
  ) %>%
  mutate(
    strategy = factor(
      strategy,
      levels = c("dual", "list", "agar"),
      labels = c(
        "Dual strategy",
        "Edible list based",
        "Agaricomycetes based"
      )
    )
  )



heat_depth <- ggplot(df_depth_prj, aes(x = strategy, y = cohort, fill = percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", percent)), size = 3) +
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    name = "%"
  ) +
  scale_x_discrete(expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(mult = 0)) +
  theme_bw() +
  labs(
    x = "Filtering strategy of marcofungi",
    y = "Sub-project") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1))


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
ggsave(file.path(path, "heat_macrofungi_depth.pdf"), heat_depth, width = 6, height = 12)




df_obs_prj <- df %>%
  mutate(
    dual = obs_macro_dual / obs_all * 100,
    list = obs_macro_list / obs_all * 100,
    agar = obs_macro_agar / obs_all * 100
  ) %>%
  group_by(cohort) %>%
  summarise(
    dual = mean(dual, na.rm = TRUE),
    list = mean(list, na.rm = TRUE),
    agar = mean(agar, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(dual, list, agar),
    names_to = "strategy",
    values_to = "percent"
  ) %>%
  mutate(
    strategy = factor(
      strategy,
      levels = c("dual", "list", "agar"),
      labels = c(
        "Dual strategy",
        "Edible list based",
        "Agaricomycetes based"
      )
    )
  )


heat_obs <- ggplot(df_obs_prj, aes(x = strategy, y = cohort, fill = percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", percent)), size = 3) +
  scale_fill_gradient(
    low = "white",
    high = "firebrick",
    name = "%"
  ) +
  scale_x_discrete(expand = expansion(mult = 0)) +
  scale_y_discrete(expand = expansion(mult = 0)) +
  theme_bw() +
  labs(
    x = "Filtering strategy of marcofungi",
    y = "Sub-project") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1))


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
ggsave(file.path(path, "heat_macrofungi_obs.pdf"), heat_obs, width = 6, height = 12)







ps_all <- readRDS("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check/ps_all.rds")
ps_fungi <- subset_taxa(ps_all, Kingdom == "k__Fungi")




mushroom_list <- read.csv("/data/wangxinyu/Tools_Dev/MycoGAP/mycogap/ref/mushroom_genus.csv")
mushroom_list <- mushroom_list$Genus
mushroom_list <- paste0("g__", mushroom_list)
head(mushroom_list)

# Find ASVs whose genus name matches entries in the mushroom database
taxa_table <- as.data.frame(ps_fungi@tax_table)
taxa_table_mushroom1 <- taxa_table %>% filter(Genus %in% mushroom_list)
head(taxa_table_mushroom1)
nrow(taxa_table_mushroom1) # 13077


# Find ASVs in taxa_table where class is Agaricomycetes
taxa_table_mushroom2 <- taxa_table %>% filter(Class == "c__Agaricomycetes")
head(taxa_table_mushroom2)
nrow(taxa_table_mushroom2) # 16816

# Remove identified ASVs
ASV_to_drop <- union(rownames(taxa_table_mushroom1), rownames(taxa_table_mushroom2))
head(ASV_to_drop)
length(ASV_to_drop) # 18988


ps_microfungi_dual <- prune_taxa(!(taxa_names(ps_fungi) %in% ASV_to_drop), ps_fungi)
ps_microfungi_list <- prune_taxa(!(taxa_names(ps_fungi) %in% rownames(taxa_table_mushroom1)), ps_fungi)
ps_microfungi_agar <- prune_taxa(!(taxa_names(ps_fungi) %in% rownames(taxa_table_mushroom2)), ps_fungi)


ps_microfungi_dual_fdp <- prune_samples(sample_sums(ps_microfungi_dual) >= 10000, ps_microfungi_dual)



samples_keep <- sample_names(ps_microfungi_dual_fdp)
length(samples_keep) # 37417
ps_microfungi_list_fdp <- prune_samples(samples_keep, ps_microfungi_list)
ps_microfungi_agar_fdp <- prune_samples(samples_keep, ps_microfungi_agar)

rm(ps_microfungi_dual, ps_microfungi_list, ps_microfungi_agar)


ps_microfungi_dual_fdp <- filter_taxa(ps_microfungi_dual_fdp, function(otu) sum(otu) > 1, TRUE)
ps_microfungi_list_fdp <- filter_taxa(ps_microfungi_list_fdp, function(otu) sum(otu) > 1, TRUE)
ps_microfungi_agar_fdp <- filter_taxa(ps_microfungi_agar_fdp, function(otu) sum(otu) > 1, TRUE)



path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
saveRDS(ps_microfungi_dual_fdp, file = file.path(path, "ps_microfungi_dual_fdp.rds"))
saveRDS(ps_microfungi_list_fdp, file = file.path(path, "ps_microfungi_list_fdp.rds"))
saveRDS(ps_microfungi_agar_fdp, file = file.path(path, "ps_microfungi_agar_fdp.rds"))



div_dual <- estimate_richness(ps_microfungi_dual_fdp, split = TRUE, measures = c("Observed", "Shannon", "Simpson")) %>% rename_with(~ paste0(.x, "_micro_dual"))
div_list <- estimate_richness(ps_microfungi_list_fdp, split = TRUE, measures = c("Observed", "Shannon", "Simpson")) %>% rename_with(~ paste0(.x, "_micro_list"))
div_agar <- estimate_richness(ps_microfungi_agar_fdp, split = TRUE, measures = c("Observed", "Shannon", "Simpson")) %>% rename_with(~ paste0(.x, "_micro_agar"))
head(div_list)




res_list <- list(
  div_dual,
  div_list,
  div_agar
)

div <- res_list %>%
  map(~ rownames_to_column(.x, "ID")) %>%
  reduce(full_join, by = "ID")

head(div)
dim(div)


# metadata
div_all <- read.csv("/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check/div_all.csv", row.names = 1)
div_all <- div_all %>% rename_with(~ paste0(.x, "_all")) %>% rownames_to_column("ID")

meta <- as.matrix(sample_data(ps_all)) %>%
  as.data.frame() %>%
  select(cohort) %>%
  rownames_to_column("ID")


meta_all <- merge(meta, div_all, by = "ID")
head(meta_all)


df <- merge(meta_all, div, by = "ID", all.y = T)
head(df)
nrow(df) # 37417

path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
write.csv(df, file.path(path, "macrofungi_filter_div_sample.csv"))



df_long <- df %>%
  pivot_longer(
    cols = -c(ID, cohort),
    names_to = c("Index", "Strategy"),
    names_pattern = "(Observed|Shannon|Simpson)_(all|micro_dual|micro_list|micro_agar)",
    values_to = "Value"
  ) %>%
  mutate(
    Strategy = factor(
      Strategy,
      levels = c("all", "micro_dual", "micro_list", "micro_agar"),
      labels = c(
        "All species",
        "Microfungi (dual)",
        "Microfungi (list)",
        "Microfungi (agar)"
      )
    ),
    Index = recode(
      Index,
      Observed = "Observed features (log10)",
      Shannon = "Shannon",
      Simpson = "Simpson"
    ),
    Index = factor(
      Index,
      levels = c("Observed features (log10)", "Shannon", "Simpson")
    ),
    Value_plot = ifelse(Index == "Observed features (log10)", log10(Value), Value)
  )

my_comparisons <- list(
  c("All species", "Microfungi (dual)"),
  c("All species", "Microfungi (list)"),
  c("All species", "Microfungi (agar)")
)

my_colors <- c(
  "All species" = "#ad6a6c",
  "Microfungi (dual)" = "#1e1e24",
  "Microfungi (list)" = "#3a3a44",
  "Microfungi (agar)" = "#565666"
)

box_sample <- ggplot(
  df_long,
  aes(x = Strategy, y = Value_plot, fill = Strategy)
) +
  geom_boxplot(
    alpha = 0.8,
    color = "black",
    outlier.shape = 16,
    outlier.size = 1
  ) +
  facet_wrap(~ Index, scales = "free_y") +
  scale_fill_manual(values = my_colors) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(after_stat(y), 2)),
      angle = 90,
    vjust = -0.6,
    size = 3,
    color = "black"
  ) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = my_comparisons,
    label = "p.format",
    size = 3,
    tip.length = 0,
    bracket.size = 0.5,
    step.increase = 0.18
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #axis.title.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Sample-based comparison",
    title = NULL
  )



path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
ggsave(file.path(path, "box_macrfungi_alpha_sample.pdf"), box_sample, width = 6, height = 3)



# PRJ based

df_project <- df %>%
  group_by(cohort) %>%
  summarise(
    across(
      c(
        Observed_all,
        Shannon_all,
        Simpson_all,
        Observed_micro_dual,
        Shannon_micro_dual,
        Simpson_micro_dual,
        Observed_micro_list,
        Shannon_micro_list,
        Simpson_micro_list,
        Observed_micro_agar,
        Shannon_micro_agar,
        Simpson_micro_agar
      ),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

head(df_project)

path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
write.csv(df_project, file.path(path, "macrofungi_filter_div_project.csv"))



df_project_long <- df_project %>%
  pivot_longer(
    cols = -cohort,
    names_to = c("Index", "Strategy"),
    names_pattern = "(Observed|Shannon|Simpson)_(all|micro_dual|micro_list|micro_agar)",
    values_to = "Value"
  ) %>%
  mutate(
    Strategy = factor(
      Strategy,
      levels = c("all", "micro_dual", "micro_list", "micro_agar"),
      labels = c(
        "All species",
        "Microfungi (dual)",
        "Microfungi (list)",
        "Microfungi (agar)"
      )
    ),
    Index = recode(
      Index,
      Observed = "Observed features (log10)",
      Shannon = "Shannon",
      Simpson = "Simpson"
    ),
    Index = factor(
      Index,
      levels = c("Observed features (log10)", "Shannon", "Simpson")
    ),
    Value_plot = ifelse(Index == "Observed features (log10)", log10(Value), Value)
  )

my_comparisons <- list(
  c("All species", "Microfungi (dual)"),
  c("All species", "Microfungi (list)"),
  c("All species", "Microfungi (agar)")
)

my_colors <- c(
  "All species" = "#ad6a6c",
  "Microfungi (dual)" = "#1e1e24",
  "Microfungi (list)" = "#3a3a44",
  "Microfungi (agar)" = "#565666"
)


box_project <- ggplot(
  df_project_long,
  aes(x = Strategy, y = Value_plot, fill = Strategy)
) +
  geom_boxplot(
    alpha = 0.8,
    color = "black",
    outlier.shape = 16,
    outlier.size = 1
  ) +
  facet_wrap(~ Index, scales = "free_y") +
  scale_fill_manual(values = my_colors) +


  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(after_stat(y), 2)),
      angle = 90,
    vjust = -0.6,
    size = 3,
    color = "black"
  ) +
  stat_compare_means(
    method = "t.test",
    paired = TRUE,
    comparisons = my_comparisons,
    label = "p.format",
    size = 3,
    tip.length = 0,
    bracket.size = 0.5,
    step.increase = 0.18
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    #axis.title.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Sub-project-based comparison",
    title = NULL
  )


path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
ggsave(file.path(path, "box_macrfungi_alpha_project.pdf"), box_project, width = 6, height = 3)





path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
ps_microfungi_list_fdp <- readRDS(file.path(path, "ps_microfungi_list_fdp.rds"))
ps_microfungi_agar_fdp <- readRDS(file.path(path, "ps_microfungi_agar_fdp.rds"))


ps_micro_list_gen <- aggregate_rare(ps_microfungi_list_fdp, level = "Genus", detection = 0, prevalence = 0.001)
ps_micro_agar_gen <- aggregate_rare(ps_microfungi_agar_fdp, level = "Genus", detection = 0, prevalence = 0.001)
ps_micro_list_gen # 741 taxa and 37417 samples
ps_micro_agar_gen # 713 taxa and 37417 samples

otu_micro_list_gen <- otu_table(ps_micro_list_gen)
otu_micro_agar_gen <- otu_table(ps_micro_agar_gen)

otu_micro_list_gen <- t(otu_micro_list_gen)
otu_micro_agar_gen <- t(otu_micro_agar_gen)


dist_micro_list_gen <- vegdist(otu_micro_list_gen, MARGIN = 1, method = "robust.aitchison")
dist_micro_agar_gen <- vegdist(otu_micro_agar_gen, MARGIN = 1, method = "robust.aitchison")
str(dist_micro_list_gen)
str(dist_micro_agar_gen)



path <- "/data/wangxinyu/Fungi_Atlas/Analysis/Revision/macrofungi_check"
saveRDS(dist_micro_list_gen, file.path(path, "dist_micro_list_gen.rds"))
saveRDS(dist_micro_agar_gen, file.path(path, "dist_micro_agar_gen.rds"))
