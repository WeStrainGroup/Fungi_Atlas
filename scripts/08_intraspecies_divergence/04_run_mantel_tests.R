#!/usr/bin/env Rscript
# Purpose: Convert SH trees to sample genetic-distance matrices and test correlation with host geographic distance.

options(echo = TRUE)

library(ape)
library(vegan)
library(tidyverse)
library(data.table)

# Read a large OTU table and restore its first column as row names.
read_otu <- function(filepath) {
  otu <- fread(filepath, header = TRUE)
  otu <- as.data.frame(otu)
  rownames(otu) <- otu[[1]]
  otu <- otu[, -1]
  return(otu)
}

# Function: compute Mantel correlation for a given SH
mantel_for_SH <- function(SH_name, otu, cluster_index, tree_dir, meta, dis_region, permutations = 999, parallel = 1) {
  cat("\nRunning:", SH_name)
  # 1. Retrieve all ASVs that belong to this SH
  SH_asvs <- cluster_index %>% filter(SH == SH_name) %>% pull(ASV) %>% as.character()

  # 2. Subset OTU to these ASVs and drop samples with zero total counts
  otu_subset <- otu %>% select(all_of(SH_asvs)) %>% filter(rowSums(across(everything())) != 0)

  # 3. Summarize per-sample dominant ASV and related counts
  otu_summary <- otu_subset %>%
    mutate(
      ID = rownames(.),
      n_ASV = rowSums(. > 0),
      dominant_ASV_count = apply(., 1, max),
      dominant_ASV = apply(., 1, function(x) colnames(otu_subset)[which.max(x)]),
      depth = rowSums(.),
      dominant_ASV_count_ratio = dominant_ASV_count / depth
    ) %>%
    select(ID, n_ASV, dominant_ASV, dominant_ASV_count, dominant_ASV_count_ratio)

  # 4. Read the SH tree and compute the ASV cophenetic distance matrix
  tree_path <- file.path(tree_dir, SH_name, paste0(SH_name, "_align.fasta.treefile"))
  tree <- read.tree(tree_path)
  tree_dis <- cophenetic.phylo(tree)

  # 5. Build mapping from sample ID to dominant ASV and initialize the distance matrix
  dominant_ASV <- unique(otu_summary$dominant_ASV)
  ID_to_ASV <- setNames(otu_summary$dominant_ASV, otu_summary$ID)
  sample_ids <- names(ID_to_ASV)

  # 6. Fill the phylogenetic distance matrix using vectorized indexing
  ASV_per_sample <- ID_to_ASV[sample_ids] # Fetch dominant ASV for each sample in the order of sample_ids
  id_dis_phylo <- tree_dis[ASV_per_sample, ASV_per_sample, drop = FALSE]

  rownames(id_dis_phylo) <- sample_ids
  colnames(id_dis_phylo) <- sample_ids

  # 7. Prepare the geographic distance matrix
  geo_info <- meta %>% select(ID, Region) %>% filter(ID %in% sample_ids)
  id_to_region <- setNames(geo_info$Region, geo_info$ID)
  region_per_sample <- id_to_region[sample_ids]  # Ordered by sample_ids

  # Build region-to-region distance matrix
  region_dist_mat <- reshape2::acast(dis_region, Region1 ~ Region2, value.var = "Distance_km")

  # Check which regions exist in the region distance matrix
  valid_id <- region_per_sample %in% rownames(region_dist_mat)

  # Filter out samples whose regions are absent from the distance matrix
  sample_ids_valid <- sample_ids[valid_id]
  region_per_sample_valid <- region_per_sample[valid_id]

  # Construct the sample-by-sample geographic distance matrix
  id_dis_geogra <- region_dist_mat[region_per_sample_valid, region_per_sample_valid, drop = FALSE]
  rownames(id_dis_geogra) <- sample_ids_valid
  colnames(id_dis_geogra) <- sample_ids_valid


  # 9. Align the two matrices (intersect sample IDs and ensure identical order)
  common_ids <- intersect(rownames(id_dis_phylo), rownames(id_dis_geogra))
  id_dis_phylo_aligned <- id_dis_phylo[common_ids, common_ids]
  id_dis_geogra_aligned <- id_dis_geogra[common_ids, common_ids]

  # Assert that row/column orders are identical
  stopifnot(identical(rownames(id_dis_phylo_aligned), rownames(id_dis_geogra_aligned)))
  stopifnot(identical(colnames(id_dis_phylo_aligned), colnames(id_dis_geogra_aligned)))

  # 10. Run the Mantel test
  cat("\nMantel test for:", SH_name)
  set.seed(6311)
  mantel_res <- mantel(as.dist(id_dis_phylo_aligned), as.dist(id_dis_geogra_aligned),
                       method = "pearson", permutations = permutations, parallel = parallel)

  # 12. Return a tidy data.frame of results
  data.frame(
    SH = SH_name,
    n_sample = length(common_ids),
    Rho = mantel_res$statistic,
    P = mantel_res$signif,
    Permutations = mantel_res$permutations
  )
}


# Load input data
otu <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_otu_fdp.csv")
cluster_index <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/cluster_index.csv", row.names = 1)
meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
dis_region <- read.delim("/data/huangkailang/project/fungi_ITS/data_for_analysis/distance_association/geographic_distance.tsv")
tree_dir <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/tree"

SHs_selected <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/SHs_selected.csv", row.names = 1)
SHs_selected <- SHs_selected %>% pull(SH) %>% as.character()
SHs_selected

# Iterate over SHs
results_list <- list()
for (SH in SHs_selected) {
  tryCatch({
     results_list[[SH]] <- mantel_for_SH(SH_name = SH,
                                        otu = otu,
                                        cluster_index = cluster_index,
                                        tree_dir = tree_dir,
                                        meta = meta,
                                        dis_region = dis_region,
                                        permutations = 999,
                                        parallel = 10)
    }, error = function(e) {
    message("Error in: ", SH, " - ", e$message)
    results_list[[SH]] <- NULL
  })
}


final_results <- bind_rows(results_list)
final_results

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered"
write.csv(final_results, file.path(path, "mantel_result.csv"))
