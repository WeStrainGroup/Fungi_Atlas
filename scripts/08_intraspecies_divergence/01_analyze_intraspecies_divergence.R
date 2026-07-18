# Purpose: Build global SH inputs, select diverse SHs, integrate tree/Mantel results, and generate Figure 5 statistics and plots.
# Panels: Figure 5B-E.

library(phyloseq)
library(microbiome)
library(data.table)
library(Biostrings)
library(tidyverse)
library(ape)
library(ggtree)
library(paletteer)
library(vegan)
library(future.apply)




# --- Merge raw OTU and taxa tables before distance clustering/filtering (filter during merge to avoid OOM) ---
# Function to read large OTU table
read_otu <- function(filepath) {
  otu <- fread(filepath, header = TRUE)
  otu <- as.data.frame(otu)
  rownames(otu) <- otu[[1]]
  otu <- otu[, -1]
  return(otu)
}

project_root1 <- "/data/wangxinyu/ITS_Public/Project/Public/Annotation/SE"
project_root2 <- "/data/wangxinyu/ITS_Public/Project/Public/Annotation/PE"
project_root3 <- "/data/wangxinyu/ITS_Public/Project/WeGut"

project1 <- list.dirs(project_root1, full.names = TRUE, recursive = FALSE)
project2 <- list.dirs(project_root2, full.names = TRUE, recursive = FALSE)
project3 <- list.dirs(project_root3, full.names = TRUE, recursive = FALSE)

projects <- c(project1, project2, project3)
projects
length(projects) # 83


# Start merging
ps_merged <- NULL  # Store merged result
not_merged_projects <- character()

for (proj in projects) {
  project_name <- basename(proj)
  cat("Processing:", project_name, "\n")

  seqtab_file <- file.path(proj, "1_dada2/result", paste0(project_name, "_seqtab.csv"))
  taxa_file   <- file.path(proj, "1_dada2/result", paste0(project_name, "_taxa.csv"))

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

    ps <- phyloseq(otu_table(as.matrix(seqtab), taxa_are_rows = FALSE),
                   tax_table(as.matrix(taxa)),
                   sample_data(samdf))

    # Basic filtering
    # Remove ASVs with unclear phylum annotation
    ps_filter <- subset_taxa(ps, Phylum != "NA" & Kingdom != "k__Eukaryota_kgd_Incertae_sedis")

    # Set ASV counts to 0 in samples where relative abundance is below filter_abundance
    otu <- ps_filter@otu_table
    otu[otu / rowSums(otu) < 1/10000] <- 0
    otu_table(ps_filter) <- otu # Return otu table to ps
    ps_filter <- filter_taxa(ps_filter, function(otu) sum(otu) > 1, TRUE) # Remove ASVs that are all zeros after filtering
    ps_filter <- prune_samples(sample_sums(ps_filter) > 0, ps_filter) # Remove samples without any reads

    # Fungal community related
    ps_fungi <- subset_taxa(ps_filter, Kingdom == "k__Fungi")

    # filter mushroom taxa from ps subject
    # Identify mushroom taxa
    mushroom_list <- read.csv(file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library/mycogap/ref/mushroom_genus.csv"))
    mushroom_list <- mushroom_list$Genus
    mushroom_list <- paste0("g__", mushroom_list)

    # Find ASVs matching genus names in the reference list
    taxa_table <- as.data.frame(ps_fungi@tax_table)
    taxa_table_mushroom1 <- taxa_table %>% filter(Genus %in% mushroom_list)

    # Find ASVs with class = Agaricomycetes in taxa_table
    taxa_table_mushroom2 <- taxa_table %>% filter(Class == "c__Agaricomycetes")

    # Remove identified ASVs
    ASV_to_drop <- union(rownames(taxa_table_mushroom1), rownames(taxa_table_mushroom2))
    ps_fungi2 <- prune_taxa(!(taxa_names(ps_fungi) %in% ASV_to_drop), ps_fungi)

    # Remove samples without any reads
    ps_fungi2 <- prune_samples(sample_sums(ps_fungi2) > 0, ps_fungi2)

    # add ASV tags
    ps_fungi2 <- add_refseq(ps_fungi2, tag = paste0(project_name, "_ASV"))

    if (is.null(ps_merged)) {
      ps_merged <- ps_fungi2
    } else {
      ps_merged <- merge_phyloseq(ps_merged, ps_fungi2)
      rm(seqtab, taxa, samdf, ps, ps_filter, ps_fungi, ps_fungi2)
      gc()
    }
  } else {
    cat("No file found: ", project_name, "\n")
    not_merged_projects <- c(not_merged_projects, project_name)
  }
}


# check
ps_merged
not_merged_projects
otu_table(ps_merged)[1:5, 1:5]
tax_table(ps_merged)[1:5, 1:5]


# Remove samples with insufficient reads
ps_merged2 <- prune_samples(sample_sums(ps_merged) >= 10000, ps_merged)
# Remove ASVs that are all zeros after filtering
ps_merged2 <- filter_taxa(ps_merged2, function(otu) sum(otu) > 1, TRUE)
ps_merged2


# Write out
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance"
project <- "unc_"

saveRDS(ps_merged, file = file.path(path, "ps_unc.rds"))
saveRDS(ps_merged2, file = file.path(path, "ps_unc_fdp.rds"))
ps_merged2 <- readRDS("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/ps_unc_fdp.rds")
ps_merged2
write.csv(readcount(ps_merged2), file.path(path, paste0(project, "readcount_fdp.csv")))
writeXStringSet(ps_merged2@refseq, file = file.path(path, paste0(project, "refseq_fdp.fasta")), format = "fasta")
write.csv(ps_merged2@tax_table, file.path(path, paste0(project, "taxa_fdp.csv")))
write.csv(ps_merged2@otu_table, file.path(path, paste0(project,"otu_fdp.csv")))






# --- Start clustering ASVs: split refseq into ITS1 and ITS2 and cluster separately ---
ASVs <- readDNAStringSet("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_refseq_fdp.fasta",
                         format = "fasta")

ASVs_ITS1 <- ASVs[grepl("ITS1", names(ASVs), ignore.case = F)]
ASVs_ITS2 <- ASVs[grepl("ITS2", names(ASVs), ignore.case = F)]

ASVs
ASVs_ITS1
ASVs_ITS2

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw"
writeXStringSet(ASVs_ITS1, file = file.path(path, "unc_refseq_fdp_ITS1.fasta"), format = "fasta")
writeXStringSet(ASVs_ITS2, file = file.path(path, "unc_refseq_fdp_ITS2.fasta"), format = "fasta")


# --- run vsearch ---
# Run 02_cluster_global_shs.sh, then resume this script.

# --- Organize and summarize clusters ---
# ITS1
cluster_index <- read.delim("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS1_v98.5.uc", header = F)
head(cluster_index)

cluster_index1 <- cluster_index %>%
        filter(V1 != "C") %>%
        mutate(V10 = ifelse(V10 == "*", V9, V10)) %>%
        select(V9, V10) %>%
        dplyr::rename(ASV = V9, rASV = V10) %>%
        add_count(rASV, name = "n_ASV") %>%
        arrange(-n_ASV) %>%
        mutate(marker = "ITS1")
head(cluster_index1)

nrow(cluster_index1)
length(unique(cluster_index1$rASV))

# ITS2
cluster_index <- read.delim("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS2_v98.5.uc", header = F)
head(cluster_index)

cluster_index2 <- cluster_index %>%
        filter(V1 != "C") %>%
        mutate(V10 = ifelse(V10 == "*", V9, V10)) %>%
        select(V9, V10) %>%
        dplyr::rename(ASV = V9, rASV = V10) %>%
        add_count(rASV, name = "n_ASV") %>%
        arrange(-n_ASV) %>%
        mutate(marker = "ITS2")
head(cluster_index2)

nrow(cluster_index2)
length(unique(cluster_index2$rASV))

cluster_index_merged <- rbind(cluster_index1, cluster_index2)
head(cluster_index_merged)
tail(cluster_index_merged)

# Assign nickname/ID to each cluster: SH1, SH2, ...
cluster_index_merged <- cluster_index_merged %>%
  group_by(marker, rASV) %>%
  mutate(SH = paste0(marker, "_SH", cur_group_id())) %>%
  ungroup()
head(cluster_index_merged)
tail(cluster_index_merged)
length(unique(cluster_index_merged$SH)) #39272

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch"
write.csv(cluster_index_merged, file.path(path, "cluster_index.csv"))




# --- Extract all ASVs for each SH into a FASTA file ---
cluster_index <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/cluster_index.csv",
                          row.names = 1)
head(cluster_index)

refseq <- readDNAStringSet("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_refseq_fdp.fasta", format = "fasta")
refseq

# Get all SHs and their ASVs
sh_to_asv_list <- cluster_index %>%
  group_by(SH) %>%
  summarise(ASVs = list(ASV), .groups = "drop")
sh_to_asv_list

# Extract sequences for each SH and save as FASTA
output_dir <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/clusters"

for (i in seq_len(nrow(sh_to_asv_list))) {
  sh_name <- sh_to_asv_list$SH[i]
  asvs <- sh_to_asv_list$ASVs[[i]]
  seqs <- refseq[names(refseq) %in% asvs]
  writeXStringSet(seqs, file = file.path(output_dir, paste0(sh_name, ".fasta")))
}





# --- Merge ASVs in the ps object based on vsearch results for cluster-level stats ---
otu <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_otu_fdp.csv")
otu <- as.data.frame(t(otu))
otu <- otu %>% rownames_to_column("ASV")
otu[1:5, 1:5]
nrow(otu)


# Add index information to seqtab (following MycoGAP code)
seqtab7 <- merge(cluster_index_merged, otu, by = "ASV", all.y = TRUE)
seqtab7 <- seqtab7 %>%
  select(-ASV, -n_ASV, -marker)

# Check
head(colnames(seqtab7))
head(seqtab7$rASV)

# Merge identical rASVs
setDT(seqtab7)

seqtab8 <- seqtab7[, lapply(.SD, sum, na.rm = TRUE), by = rASV, .SDcols = where(is.numeric)]
seqtab8 <- as.data.frame(seqtab8)

# Convert sequence to row names and transpose
seqtab8 <- seqtab8 %>% column_to_rownames(var = "rASV")

seqtab8[1:3, 1:3]

seqtab_cluster <- t(seqtab8)
seqtab_cluster[1:3, 1:3]
dim(seqtab_cluster)


taxa <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_taxa_fdp.csv", row.names = 1)
head(taxa)
nrow(taxa)

taxa_cluster <- taxa %>%
  rownames_to_column("ASV") %>%
  filter(ASV %in% colnames(seqtab_cluster)) %>%
  column_to_rownames("ASV")

head(taxa_cluster)
dim(taxa_cluster)



# Rename each cluster to SH
# Build a name map
name_map <- cluster_index_merged %>%
  distinct(rASV, SH) %>%
  deframe()  # convert to named vector
head(name_map)

# otu
seqtab_cluster <- seqtab_cluster %>%
  as.data.frame() %>%
  rename_with(~ ifelse(.x %in% names(name_map), name_map[.x], .x))

seqtab_cluster[1:3, 1:3]


taxa_cluster <- taxa_cluster %>%
  rownames_to_column("rASV") %>%
  mutate(rASV = ifelse(rASV %in% names(name_map), name_map[rASV], rASV)) %>%
  column_to_rownames("rASV")

head(taxa_cluster)


# Write out
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered"
write.csv(seqtab_cluster, file.path(path, "clu_otu_fdp.csv"))
write.csv(taxa_cluster, file.path(path, "clu_taxa_fdp.csv"))



# Build a ps object
# Construct the phyloseq subject
otu_table <- phyloseq::otu_table(as.matrix(seqtab_cluster), taxa_are_rows = FALSE)

sam_data <- data.frame(row.names = rownames(otu_table), ID = rownames(otu_table))

meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
sam_data <- merge(sam_data, meta, by = "ID", all.x = T)
rownames(sam_data) <- sam_data$ID
sam_data <- sample_data(sam_data)
head(sam_data)

tax_table <- phyloseq::tax_table(as.matrix(tax_cluster))

refseq1 <- readDNAStringSet("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS1_v98.5.fasta", format = "fasta")
refseq2 <- readDNAStringSet("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS2_v98.5.fasta", format = "fasta")
refseq <- c(refseq1, refseq2)
# Rename sequences
  setdiff(names(refseq), names(name_map)) # check for inconsistencies
  names(refseq) <- name_map[names(refseq)]

ps <- phyloseq(otu_table,
               sam_data,
               tax_table,
               refseq)
ps

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered"
saveRDS(ps, file = file.path(path, "ps_clu_fdp.rds"))





# --- Compute prevalence and add to SH metadata ---
ps_ITS1 <- subset_samples(ps, PCR_marker == "ITS1")
ps_ITS2 <- subset_samples(ps, PCR_marker == "ITS2")

ps_ITS1
ps_ITS2

# Remove ASVs that are all zeros after filtering
ps_ITS1 <- filter_taxa(ps_ITS1, function(otu) sum(otu) > 0, TRUE)
ps_ITS2 <- filter_taxa(ps_ITS2, function(otu) sum(otu) > 0, TRUE)

ps_ITS1
ps_ITS2


prev_ITS1 <- as.data.frame(prevalence(ps_ITS1, sort = TRUE, detection = 0))
colnames(prev_ITS1) <- "prev"
prev_ITS1 <- prev_ITS1 %>% rownames_to_column("SH")
head(prev_ITS1)
nrow(prev_ITS1) #9289

prev_ITS2 <- as.data.frame(prevalence(ps_ITS2, sort = TRUE, detection = 0))
colnames(prev_ITS2) <- "prev"
prev_ITS2 <- prev_ITS2 %>% rownames_to_column("SH")
head(prev_ITS2)
nrow(prev_ITS2) #29966

# Merge
prev <- rbind(prev_ITS1, prev_ITS2)
nrow(prev) # 39255

# Merge with cluster_index_merged
head(cluster_index_merged)
cluster_map <- cluster_index_merged %>%
    select(-ASV) %>%
    distinct() %>%
    as.data.frame()
head(cluster_map)

cluster_map2 <- merge(cluster_map, prev, by = "SH", all = T)

cluster_map2 <- cluster_map2 %>%
    arrange(-n_ASV)

head(cluster_map2)
nrow(cluster_map2) # 39272

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered"
write.csv(cluster_map2, file.path(path, "meta_SH.csv"))




# --- Visualize SH distribution ---
meta_SH <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/meta_SH.csv",
       row.names =1)

head(meta_SH)

# Check NA (likely samples missing in metadata)
table(is.na(meta_SH$prev))
na_rows <- meta_SH %>% filter(is.na(prev)) # All have n_ASV = 1; can ignore for plotting
na_rows

threshold <- 250
SHs_include <- meta_SH %>%
  filter(n_ASV >= threshold)

n_SHs_include <- nrow(SHs_include)
n_SHs_include # 67

n_SHs_notsingle <- meta_SH %>%
  filter(n_ASV > 1) %>%
  nrow()
n_SHs_notsingle #17716 for 1; 39272 for all

col <- c("ITS1" = "#bc4749",
         "ITS2" = "#457b9d")

meta_SH_shuffled <- meta_SH %>%
  slice_sample(prop = 1)  # randomly shuffle all rows

p <- ggplot(meta_SH_shuffled, aes(x = prev, y = n_ASV, color = marker)) +
        geom_point(alpha = 0.5, shape = 16, size = 2.5) +
        labs(x = "Prevalence of SHs",
             y = "# ASVs clustered into each SHs",
             color = NULL) +
        scale_y_log10(limits = c(1, 10000)) +
        xlim(0,1) +
        scale_color_manual(values = col) +
        geom_hline(yintercept = threshold, linewidth = 0.5, linetype = "dashed", color = "black") +
        annotate("text", x = 0.6, y = threshold - 100,
                 label = paste0("threshold: ", 250), size = 3, hjust = 0) +
        annotate("text", x = 0, y = 5000,
                 label = paste0(n_SHs_include, " SHs were included"),
                 size = 3, hjust = 0) +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", linewidth = 1),
              legend.background = element_rect(fill = "transparent", color = NA),
              legend.key = element_rect(fill = "transparent", color = NA),
              legend.position = c(0.8, 0.1),
              legend.key.size = unit(0.15, "in"))

# No ITS1/ITS2 distinction
p <- ggplot(meta_SH_shuffled, aes(x = prev, y = n_ASV)) +
        geom_point(alpha = 0.2, shape = 16, size = 2.5, color = "black") +
        labs(x = "Prevalence of SHs",
             y = "# ASVs clustered into each SHs",
             color = NULL) +
        scale_y_log10(limits = c(1, 10000)) +
        xlim(0,1) +
        geom_hline(yintercept = threshold, linewidth = 0.5, linetype = "dashed", color = "black") +
        annotate("text", x = 0.5, y = threshold - 100,
                 label = paste0("threshold: ", 250), hjust = 0) +
        annotate("text", x = 0, y = 5000,
                 label = paste0(n_SHs_include, " SHs were included"),
                 hjust = 0) +
        theme_bw() +
        theme(panel.border = element_rect(color = "black", linewidth = 1),
              legend.background = element_rect(fill = "transparent", color = NA),
              legend.key = element_rect(fill = "transparent", color = NA),
              legend.position = c(0.8, 0.1),
              legend.key.size = unit(0.15, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/figure"
ggsave(file.path(path, "scatter_SH_meta.pdf"), p, width = 3, height = 3)

write.csv(SHs_include, file.path(path, "SHs_selected.csv"))


# --- Run 03_build_sh_phylogenies.sh to build phylogenetic trees for selected SHs ---

# --- Run 04_run_mantel_tests.R to prepare data and perform Mantel tests, then resume this script ---


# --- Tidy Mantel test results ---
res_mantel <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/mantel_result.csv", row.names = 1)
head(res_mantel)

SHs_selected <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/SHs_selected.csv", row.names = 1)
head(SHs_selected)

taxa <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/clu_taxa_fdp.csv", row.names = 1)
taxa <- taxa %>% rownames_to_column("SH")
head(taxa)

# Merge metadata
res_mantel2 <- merge(res_mantel, SHs_selected, by = "SH")
head(res_mantel2)

res_mantel3 <- merge(res_mantel2, taxa, by = "SH", all = F)
head(res_mantel3)
nrow(res_mantel3)

# Arrange order/columns
res_mantel4 <- res_mantel3 %>%
   mutate(P_adj = p.adjust(P, method = "fdr")) %>%
   select(SH, marker, n_sample, Rho, P, P_adj, Permutations, prev, n_ASV, rASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
   arrange(P_adj, -Rho)

head(res_mantel4)


# Write out
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered"
write.csv(res_mantel4, file.path(path, "mantel_result_final.csv"))




# --- Visualize Mantel test results ---
res_mantel <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/mantel_result_final.csv", row.names = 1)
res_mantel <- res_mantel %>% filter(P_adj <= 0.05) %>% select(SH, n_sample, Rho, P_adj, prev, n_ASV, Genus, Species)
res_mantel <- res_mantel %>%
  mutate(
    Genus = sub("^g__", "", Genus),
    Species = sub("^s__", "", Species),
    Species = ifelse(is.na(Species) | Species == "", "sp.", Species),
    taxa = paste0(Genus, " ", Species)
  ) %>%
  select(-Genus, -Species) %>%
  arrange(desc(Rho))

head(res_mantel)
table(res_mantel$n_ASV)
table(res_mantel$n_sample)


# Split labels
label_down <- res_mantel %>% slice(1:2) %>% mutate(label_y = Rho - 0.02)
label_up <- res_mantel %>% slice(-c(1,2)) %>% mutate(label_y = Rho + 0.02)


bubble <- ggplot(res_mantel, aes(y = Rho, x = reorder(SH, Rho))) +
    geom_point(aes(size = n_sample, color = n_ASV), alpha = 0.95, shape = 16) +
    geom_text(data = label_down, aes(y = label_y, label = taxa), size = 3, hjust = 1, fontface = "italic", angle = 90) +
    geom_text(data = label_up, aes(y = label_y, label = taxa), size = 3, hjust = 0, fontface = "italic", angle = 90) +
  scale_color_gradient(low = "#98c1d9", high = "#000814", name = "# ASVs") +
  scale_size(name = "# Samples",
            breaks = c(5000, 15000, 25000),
            labels = c(5000, 15000, 25000)) +
  labs(y = "Rho (geo. dist and phylo. dist)",
       x = "Species hypotheses (Mantel test; FDR < 0.05)") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.6)) +
  theme(panel.border = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA),
            legend.position = "right",
            legend.key.size = unit(0.15, "in"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/figure"
ggsave(file.path(path, "bubble_mantel.pdf"), bubble, width = 9, height = 3, unit = "in")



# --- For each sample and selected SH, compute the dominant ASV ---
otu <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_otu_fdp.csv")
cluster_index <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/cluster_index.csv", row.names = 1)

SHs_selected <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/SHs_selected.csv", row.names = 1)
SHs_selected <- SHs_selected %>% pull(SH) %>% as.character()
SHs_selected


process_SH <- function(SH_name) {
  # Get all ASVs for this SH
  SH_asvs <- cluster_index %>% filter(SH == SH_name) %>% pull(ASV) %>% as.character()

  # Subset OTU columns for this SH and drop empty samples
  otu_subset <- otu %>% select(all_of(SH_asvs)) %>% filter(rowSums(across(everything())) != 0)

  # If no samples, return NULL
  if (nrow(otu_subset) == 0) return(NULL)

  # Summarize dominant ASV per sample and counts
  otu_summary <- otu_subset %>%
    mutate(
      ID = rownames(.),
      n_ASV = rowSums(. > 0),
      dominant_ASV_count = apply(., 1, max),
      dominant_ASV = apply(., 1, function(x) colnames(otu_subset)[which.max(x)]),
      depth = rowSums(.),
      dominant_ASV_count_ratio = dominant_ASV_count / depth
    ) %>%
    select(ID, n_ASV, dominant_ASV, dominant_ASV_count, dominant_ASV_count_ratio) %>%
    mutate(SH = SH_name)

  rm(otu_subset)
  gc()
  return(otu_summary)
}

# Increase global size limit to 30GB; OTU too large for parallel otherwise
options(future.globals.maxSize = 30 * 1024^3)
future::plan(multicore, workers = 8)
result_list <- list()
result_list <- future_lapply(SHs_selected, process_SH)
length(result_list)

final_result <- bind_rows(result_list)
rownames(final_result) <- NULL
head(final_result)
nrow(final_result)


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered"
write.csv(final_result, file.path(path, "sample_dominant_ASV.csv"))


# --- Visualize per-sample ASV counts per SH and dominant_ASV_count_ratio ---
head(final_result)

# Compute medians
median_nASV <- round(median(final_result$n_ASV), 0)
median_ratioASV <- round(median(final_result$dominant_ASV_count_ratio), 0)

median_nASV # 1
median_ratioASV # 1

prop_nASV <- round(mean(final_result$n_ASV == 1)*100, 1)
prop_nASV

prop_ratioASV <-round(mean(final_result$dominant_ASV_count_ratio > 0.9)*100, 1)
prop_ratioASV


his1 <- ggplot(final_result, aes(x = n_ASV)) +
  geom_histogram( position = "identity", alpha = 0.5, color = "black", linewidth = 0.25) +
      scale_x_log10() +
      labs(x = "# ASV per SH per sample",
           y = "# Samples",
           fill = NULL) +
      annotate("text", x = 15 , y = 35e4,
               label = paste0("# ASV = 1: ", prop_nASV, "%")) +
      theme_bw() +
      theme(panel.border = element_rect(color = "black", linewidth = 1))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/figure"
ggsave(file.path(path, "his_n_ASV.pdf"), his1, width = 3, height = 3)


his2 <- ggplot(final_result, aes(x = dominant_ASV_count_ratio)) +
  geom_histogram( position = "identity", alpha = 0.5, color = "black", linewidth = 0.25) +
      labs(x = "Ratio of dominant ASV per SH",
           y = "# Samples",
           fill = NULL) +
      annotate("text", x = 0.5 , y = 35e4,
               label = paste0("ratio > 0.9: ", prop_ratioASV, "%")) +
      theme_bw() +
      theme(panel.border = element_rect(color = "black", linewidth = 1))


path <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/figure"
ggsave(file.path(path, "his_ratio_ASV.pdf"), his2, width = 3, height = 3)





# --- Phylogenetic tree visualization ---
otu_summary <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/sample_dominant_ASV.csv", row.names = 1)
meta <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_meta_final.csv", row.names = 1)
tree_dir <- "/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/tree"

SHs_selected <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/clustered/SHs_selected.csv", row.names = 1)
SHs_selected <- SHs_selected %>% pull(SH) %>% as.character()
SHs_selected

# Merge metadata
meta2 <- meta %>% select(ID, Continent, Country)
otu_summary_meta <- merge(otu_summary, meta2, by = "ID", all.x = T)
head(otu_summary_meta)
unique(otu_summary$dominant_ASV)
