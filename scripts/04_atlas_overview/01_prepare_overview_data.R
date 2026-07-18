# Purpose: Build downstream atlas tables, CLR profiles, distance matrices, rarefaction results, and cross-site ordinations.
# Panels: Figure 2A, 2C, 2D, and 2G; Figure S5A-B and S5D.

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

# --- From ASV-level merged ps: aggregate to Genus and filter to obtain analysis datasets ---
ps <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/ps_merged_ASV_filterdp.rds")

ps_gen <- aggregate_taxa(ps, level = "Genus")
ps_gen # 2294 taxa

ps_gen_0.01 <- aggregate_rare(ps, level = "Genus", detection = 0, prevalence = 0.01)
ps_gen_0.01 # 200 taxa

ps_gen_0.005 <- aggregate_rare(ps, level = "Genus", detection = 0, prevalence = 0.005)
ps_gen_0.005 # 297 taxa

ps_gen_0.001 <- aggregate_rare(ps, level = "Genus", detection = 0, prevalence = 0.001)
ps_gen_0.001 # 688 taxa

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/genus"
write.csv(otu_table(ps_gen), file = file.path(path, "otu_gen_p0_filterdp.csv"))
write.csv(tax_table(ps_gen), file = file.path(path, "tax_gen_p0_filterdp.csv"))
write.csv(otu_table(ps_gen_0.01), file = file.path(path, "otu_gen_p0.01_filterdp.csv"))
write.csv(tax_table(ps_gen_0.01), file = file.path(path, "tax_gen_p0.01_filterdp.csv"))
write.csv(otu_table(ps_gen_0.005), file = file.path(path, "otu_gen_p0.005_filterdp.csv"))
write.csv(tax_table(ps_gen_0.005), file = file.path(path, "tax_gen_p0.005_filterdp.csv"))
write.csv(otu_table(ps_gen_0.001), file = file.path(path, "otu_gen_p0.001_filterdp.csv"))
write.csv(tax_table(ps_gen_0.001), file = file.path(path, "tax_gen_p0.001_filterdp.csv"))


# --- CLR transform for Genus table with prevalence > 1% (used in main analyses) ---
otu_gen <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/genus/otu_gen_p0.01_filterdp.csv")
otu_gen[1:5, 1:5]

# CLR transform
otu_gen_clr <- vegan::decostand(otu_gen, MARGIN = 2, method = "clr", pseudocount = 1)

# Check: within each sample, the sum across taxa should be ~0
otu_gen_clr[1:5, 1:5]
colSums(otu_gen_clr[ , 1:10], na.rm = TRUE)

# Write output
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/genus"
write.csv(otu_gen_clr, file.path(path, "otu_gen_P0.01_filterdp_clr.csv"))


# --- Compute Genus-level distance matrices using prevalence > 0.1% ---
otu_gen <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/genus/otu_gen_p0.001_filterdp.csv")
otu_gen <- t(otu_gen) # taxa on the col and site on the row
otu_gen[1:5, 1:5]

# robust.aitchison
otu_gen_dist <- vegdist(otu_gen, MARGIN = 1, method = "robust.aitchison") # MARGIN = 1 computes distance between rows
str(otu_gen_dist)

# Save R object
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/distance"
saveRDS(otu_gen_dist, file.path(path, "otu_gen_0.001_filterdp_robait_dist.rds"))

# Export as readable matrix
otu_gen_dist2 <- as.matrix(otu_gen_dist)
otu_gen_dist2[1:5, 1:5]

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/distance"
write.csv(otu_gen_dist2, file.path(path, "otu_gen_0.001_filterdp_robait_dist.csv"))

# Jaccard
otu_gen_dist <- vegdist(otu_gen, MARGIN = 1, method = "jaccard") # MARGIN = 1 operates on rows
str(otu_gen_dist)

# Save R object
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/distance"
saveRDS(otu_gen_dist, file.path(path, "otu_gen_0.001_filterdp_jaccard_dist.rds"))

# Export as readable matrix
otu_gen_dist2 <- as.matrix(otu_gen_dist)
otu_gen_dist2[1:5, 1:5]

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/distance"
write.csv(otu_gen_dist2, file.path(path, "otu_gen_0.001_filterdp_jaccard_dist.csv"))


# --- From ASV-level merged ps: compute OTU and TAX tables for other taxonomic levels ---
ps <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/ps_merged_ASV_filterdp.rds")

# Build phyloseq objects at different taxonomic levels
ps_gen <- aggregate_taxa(ps, level = "Genus")
ps_fam <- aggregate_taxa(ps, level = "Family")
ps_ord <- aggregate_taxa(ps, level = "Order")
ps_cla <- aggregate_taxa(ps, level = "Class")
ps_phy <- aggregate_taxa(ps, level = "Phylum")

ps_gen # 2294 taxa
ps_fam # 682 taxa
ps_ord # 256 taxa
ps_cla # 95 taxa
ps_phy # 20 taxa

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels"
write.csv(otu_table(ps_gen), file = file.path(path, "otu_gen_filterdp.csv"))
write.csv(tax_table(ps_gen), file = file.path(path, "tax_gen_filterdp.csv"))
write.csv(otu_table(ps_fam), file = file.path(path, "otu_fam_filterdp.csv"))
write.csv(tax_table(ps_fam), file = file.path(path, "tax_fam_filterdp.csv"))
write.csv(otu_table(ps_ord), file = file.path(path, "otu_ord_filterdp.csv"))
write.csv(tax_table(ps_ord), file = file.path(path, "tax_ord_filterdp.csv"))
write.csv(otu_table(ps_cla), file = file.path(path, "otu_cla_filterdp.csv"))
write.csv(tax_table(ps_cla), file = file.path(path, "tax_cla_filterdp.csv"))
write.csv(otu_table(ps_phy), file = file.path(path, "otu_phy_filterdp.csv"))
write.csv(tax_table(ps_phy), file = file.path(path, "tax_phy_filterdp.csv"))



# --- Compute taxa prevalence at each level ----
prev_gen <- as.data.frame(prevalence(ps_gen, sort = TRUE, detection = 0))
prev_fam <- as.data.frame(prevalence(ps_fam, sort = TRUE, detection = 0))
prev_ord <- as.data.frame(prevalence(ps_ord, sort = TRUE, detection = 0))
prev_cla <- as.data.frame(prevalence(ps_cla, sort = TRUE, detection = 0))
prev_phy <- as.data.frame(prevalence(ps_phy, sort = TRUE, detection = 0))

colnames(prev_gen) <- "prevalence"
colnames(prev_fam) <- "prevalence"
colnames(prev_ord) <- "prevalence"
colnames(prev_cla) <- "prevalence"
colnames(prev_phy) <- "prevalence"

# For each level, rename 'Unknown' to prefixed 'Unassigned' to avoid merge issues
rownames(prev_gen)[rownames(prev_gen) == "Unknown"] <- "g__Unassigned"
rownames(prev_fam)[rownames(prev_fam) == "Unknown"] <- "f__Unassigned"
rownames(prev_ord)[rownames(prev_ord) == "Unknown"] <- "o__Unassigned"
rownames(prev_cla)[rownames(prev_cla) == "Unknown"] <- "c__Unassigned"
rownames(prev_phy)[rownames(prev_phy) == "Unknown"] <- "p__Unassigned"

# Add a 'level' tag
prev_gen$level <- "Genus"
prev_fam$level <- "Family"
prev_ord$level <- "Order"
prev_cla$level <- "Class"
prev_phy$level <- "Phylum"

prev_merged <- rbind(prev_phy, prev_cla, prev_ord, prev_fam, prev_gen)
table(prev_merged$level)
head(prev_merged)
tail(prev_merged)

write.csv(prev_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/prev_all_level.csv")



# --- For each level, compute the proportion of taxa with prevalence ≥ threshold ---
prev <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/prev_all_level.csv", row.names = 1)
thresholds <- seq(0, 1, by = 0.0001)
levels_list <- unique(prev$level)
levels_list

# For each level and threshold, compute the proportion of taxa with prevalence ≥ threshold
prevalence_summary <- do.call(rbind, lapply(levels_list, function(lv) {
  sub_df <- prev %>% filter(level == lv)
  data.frame(
    prevalence_threshold = thresholds,
    taxa_proportion = sapply(thresholds, function(th) {
      sum(sub_df$prevalence >= th) / nrow(sub_df)
    }),
    level = lv
  )
}))
head(prevalence_summary)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels"
write.csv(prevalence_summary, file.path(path, "prev_all_level_steps.csv"))




# --- Compute taxa abundance at each level ----
# Load OTU tables for each level
otu_gen <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_gen_filterdp.csv")
otu_fam <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_fam_filterdp.csv")
otu_ord <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_ord_filterdp.csv")
otu_cla <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_cla_filterdp.csv")
otu_phy <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_phy_filterdp.csv")

otu_gen[1:5, 1:5]

# Helper: compute total abundance across all samples (rows = taxa)
cal_abund <- function(otu, level) {
  abund <- rowSums(otu) / sum(otu)
  data.frame(row.names = rownames(otu),
             abundance = abund,
             level = level)
}

# Compute abundance by level
abund_gen <- cal_abund(otu_gen, "Genus")
abund_fam <- cal_abund(otu_fam, "Family")
abund_ord <- cal_abund(otu_ord, "Order")
abund_cla <- cal_abund(otu_cla, "Class")
abund_phy <- cal_abund(otu_phy, "Phylum")

rownames(abund_gen)[rownames(abund_gen) == "Unknown"] <- "g__Unassigned"
rownames(abund_fam)[rownames(abund_fam) == "Unknown"] <- "f__Unassigned"
rownames(abund_ord)[rownames(abund_ord) == "Unknown"] <- "o__Unassigned"
rownames(abund_cla)[rownames(abund_cla) == "Unknown"] <- "c__Unassigned"
rownames(abund_phy)[rownames(abund_phy) == "Unknown"] <- "p__Unassigned"

# Merge across levels
abund <- rbind(abund_gen, abund_fam, abund_ord, abund_cla, abund_phy)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels"
write.csv(abund, file.path(path, "abund_all_level.csv"))


# --- Rarefaction curves: number of unique taxa vs. number of samples ---
# Load OTU tables at each level
otu_phy <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_phy_filterdp.csv")
otu_cla <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_cla_filterdp.csv")
otu_ord <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_ord_filterdp.csv")
otu_fam <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_fam_filterdp.csv")
otu_gen <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/otu_gen_filterdp.csv")

otu_phy[1:5, 1:5]

# Helper: estimate unique taxa for given sample sizes
estimate_unique_taxa <- function(otu_table, step = 1000, n_rep = 1000, seed = 123) {
  set.seed(seed)
  otu_table <- as.data.frame(otu_table)
  otu_table <- otu_table[rowSums(otu_table) > 0, ]  # Remove taxa with all-zero counts
  sample_names <- colnames(otu_table)
  total_samples <- length(sample_names)

  result_list <- list()

  for (n in seq(1, total_samples, by = step)) {
    for (rep in 1:n_rep) {
      sampled <- sample(sample_names, n)
      subset <- otu_table[, sampled, drop = FALSE]
      unique_taxa <- sum(rowSums(subset) > 0)
      result_list[[length(result_list) + 1]] <- data.frame(
        sample_size = n,
        unique_taxa = unique_taxa
      )
    }
  }
  result_df <- do.call(rbind, result_list)
  return(result_df)
}

# Compute rarefaction for each level
raref_phy <- estimate_unique_taxa(otu_phy)
raref_cla <- estimate_unique_taxa(otu_cla)
raref_ord <- estimate_unique_taxa(otu_ord)
raref_fam <- estimate_unique_taxa(otu_fam)
raref_gen <- estimate_unique_taxa(otu_gen)

raref_phy$level <- "Phylum"
raref_cla$level <- "Class"
raref_ord$level <- "Order"
raref_fam$level <- "Family"
raref_gen$level <- "Genus"

raref_merged <- rbind(raref_phy, raref_cla, raref_ord, raref_fam, raref_gen)
write.csv(raref_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_levels/raref_all_level.csv")


# --- Reads-based rarefaction analysis ---
project_root1 <- "/data/wangxinyu/ITS_Public/Project/Public/Annotation/SE"
project_root2 <- "/data/wangxinyu/ITS_Public/Project/Public/Annotation/PE"
project_root3 <- "/data/wangxinyu/ITS_Public/Project/WeGut"

project1 <- list.dirs(project_root1, full.names = TRUE, recursive = FALSE)
project2 <- list.dirs(project_root2, full.names = TRUE, recursive = FALSE)
project3 <- list.dirs(project_root3, full.names = TRUE, recursive = FALSE)

projects <- c(project1, project2, project3)
projects
length(projects) # 83



# Define per-project processing function
process_proj <- function(proj) {
  project_name <- basename(proj)
  cat(project_name, "\n")
  otu_file <- file.path(proj, "3_phyloseq/fungi/ASV", paste0(project_name, "_otu_fungi_filterdp.csv"))

  if (!file.exists(otu_file)) return(NULL)

  otu <- read_otu(otu_file)
  n_sample <- min(50, nrow(otu))
  otu_sub <- otu[sample(1:nrow(otu), n_sample), ]
  otu_sub <- otu_sub[, colSums(otu_sub) > 0]
  otu_sub_t <- t(otu_sub)

  out <- tryCatch({
    out <- iNEXT(otu_sub_t, q=c(0,1,2), datatype = "abundance", knots = 50, endpoint = 20000, nboot = 100)
  }, error = function(e) return(NULL))

  if (is.null(out)) return(NULL)

  df1 <- fortify(out, type = 1)
  df2 <- fortify(out, type = 2)
  df3 <- fortify(out, type = 3)

  df1$PRJ <- project_name
  df2$PRJ <- project_name
  df3$PRJ <- project_name

  rm(otu, otu_sub, otu_sub_t, out)
  gc()

  return(list(df1 = df1, df2 = df2, df3 = df3))

}

# Parallel processing (detectCores() auto-detects threads)
result_list <- mclapply(projects, process_proj, mc.cores = 96)
head(result_list)
tail(result_list)

# Drop NULL entries
result_list <- result_list[!sapply(result_list, is.null)]

# Combine results
df1_all <- list()
df2_all <- list()
df3_all <- list()

df1_all <- do.call(rbind, lapply(result_list, `[[`, "df1"))
df2_all <- do.call(rbind, lapply(result_list, `[[`, "df2"))
df3_all <- do.call(rbind, lapply(result_list, `[[`, "df3"))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/rarefaction"
write.csv(df1_all, file.path(path, "rarefaction_reads_type1.csv"))
write.csv(df2_all, file.path(path, "rarefaction_reads_type2.csv"))
write.csv(df3_all, file.path(path, "rarefaction_reads_type3.csv"))















# --- Compute between-site distance matrices and perform MDS ---
# Distance matrices based on merged dataset (or rebuild from exported tables)
otu_gen <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged/otu_gen_p0_filterdp_merged.csv")
otu_gen[1:5, 1:5] # Rows are samples

# robust.aitchison
otu_gen_dist_roai <- vegdist(otu_gen, MARGIN = 1, method = "robust.aitchison") # MARGIN = 1 operates on rows
str(otu_gen_dist_roai)

# jaccard
otu_gen_dist_jacc <- vegdist(otu_gen, MARGIN = 1, method = "jaccard") # MARGIN = 1 operates on rows
str(otu_gen_dist_jacc)

# Output
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged"
saveRDS(otu_gen_dist_roai, file.path(path, "otu_gen_robait_dist_merged.rds"))
saveRDS(otu_gen_dist_jacc, file.path(path, "otu_gen_jaccard_dist_merged.rds"))


# Principal coordinates analysis (MDS)
otu_gen_dist <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged/otu_gen_robait_dist_merged.rds")
otu_gen_dist <- as.matrix(otu_gen_dist)

# MDS
set.seed(6311)
mds <- bigmds::divide_conquer_mds(x = otu_gen_dist, l = 300, c_points = 5 * 10, r = 10, n_cores = 16)
str(mds)

# Explained variance
eigen <- mds$eigen / sum(mds$eigen)
eigen
eigen <- as.data.frame(eigen)
rownames(eigen) <- paste0("MDS", 1:10)
head(eigen)

# Export coordinates
res <- data.frame(mds$points)
rownames(res) <- rownames(otu_gen_dist)
colnames(res) <- paste0("MDS", 1:10)
head(res)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged"
write.csv(res, file.path(path, "MDS_site_point.csv"))
write.csv(eigen, file.path(path, "MDS_site_eigen.csv"))


# Statistical tests
# --- Test whether microbiome composition differs across sites ---
otu_gen_dist <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged/otu_gen_robait_dist_merged.rds")
str(otu_gen_dist)
meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged/meta_gen_p0_filterdp_merged.csv", row.names = 1)
head(meta)

# Align sample order to match PERMANOVA requirements
meta <- meta[match(rownames(as.matrix(otu_gen_dist)), meta$id), ]
head(meta$id)
head(rownames(as.matrix(otu_gen_dist)))

# adonis2-based
adonis2 <- adonis2(otu_gen_dist ~ site , data = meta, permutations = 999, parallel = 2)
adonis2

# anosim-based
anosim <- anosim(otu_gen_dist, grouping = meta$site, permutations = 999, parallel = 2)
anosim

# Output
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged"
saveRDS(adonis2, file.path(path, "adonis2_site.rds"))
saveRDS(anosim, file.path(path, "anosim_site.rds"))
