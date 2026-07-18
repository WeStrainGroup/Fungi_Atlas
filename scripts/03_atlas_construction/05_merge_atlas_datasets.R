# Purpose: Merge project-level depth, phyloseq, diversity, metadata, mouse-gut, and environmental datasets.

#
library(phyloseq)
library(data.table)
library(Biostrings)
library(tidyverse)

# Read a large OTU table and restore its first column as row names.
read_otu <- function(filepath) {
  otu <- fread(filepath, header = TRUE)
  otu <- as.data.frame(otu)
  rownames(otu) <- otu[[1]]
  otu <- otu[, -1]
  return(otu)
}

# --- Merge depth info across projects ---
project_root1 <- "/data/wangxinyu/ITS_Public/Project/Public/Annotation/SE"
project_root2 <- "/data/wangxinyu/ITS_Public/Project/Public/Annotation/PE"
project_root3 <- "/data/wangxinyu/ITS_Public/Project/WeGut"

project1 <- list.dirs(project_root1, full.names = TRUE, recursive = FALSE)
project2 <- list.dirs(project_root2, full.names = TRUE, recursive = FALSE)
project3 <- list.dirs(project_root3, full.names = TRUE, recursive = FALSE)

projects <- c(project1, project2, project3)
projects
length(projects) # 83

depth_raw_merged <- NULL  # Container for merged results
not_merged_projects <- character()

for (proj in projects) {
  project_name <- basename(proj)
  cat("Processing:", project_name, "\n")

  depth_raw_file <- file.path(proj, "1_dada2/check", paste0(project_name, "_fq_stats_raw.tsv"))

  if (file.exists(depth_raw_file)) {
    depth_raw <- read.delim(depth_raw_file, header = T, sep = "\t")
    depth_raw$PRJ <- project_name

    if (is.null(depth_raw_merged)) {
      depth_raw_merged <- depth_raw
    } else {
      depth_raw_merged <- rbind(depth_raw_merged, depth_raw)
      rm(depth_raw)
      gc()
    }
  } else {
    cat("No file found: ", project_name, "\n")
    not_merged_projects <- c(not_merged_projects, project_name)
  }
}

not_merged_projects
head(depth_raw_merged)
tail(depth_raw_merged)

# Replace '.fastq.gz' with '.fq.gz'
depth_raw_merged$file <- sub("\\.fastq\\.gz$", ".fq.gz", depth_raw_merged$file)
head(depth_raw_merged)

# Drop rows where 'file' ends with '.2.fq.gz'
depth_raw_merged2 <- depth_raw_merged[!grepl("\\.2\\.fq\\.gz$", depth_raw_merged$file), ]
nrow(depth_raw_merged)
nrow(depth_raw_merged2)

# Remove suffixes '.fq.gz' and '.1.fq.gz' from 'file'
depth_raw_merged2$file <- sub("(\\.1)?\\.fq\\.gz$", "", depth_raw_merged2$file)
head(depth_raw_merged2)
tail(depth_raw_merged2)

# Use 'file' column as row names
rownames(depth_raw_merged2) <- depth_raw_merged2$file
head(depth_raw_merged2)
nrow(depth_raw_merged2) # 47300

# Write output
write.csv(depth_raw_merged2, "/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_raw.csv")


# Merge read counts for all, fungi, mushroom, and plant
depth_merged <- NULL
not_merged_projects <- character()

for (proj in projects) {
  project_name <- basename(proj)
  cat("Processing:", project_name, "\n")

  #depth_file <- file.path(proj, "3_phyloseq/all", paste0(project_name, "_readcount_all.csv"))
  #depth_file <- file.path(proj, "3_phyloseq/fungi/ASV", paste0(project_name, "_readcount_fungi.csv"))
  #depth_file <- file.path(proj, "3_phyloseq/vegetable/plant/ASV", paste0(project_name, "_readcount_plant.csv"))
  depth_file <- file.path(proj, "3_phyloseq/vegetable/mushroom/ASV", paste0(project_name, "_readcount_mushroom.csv"))


  if (file.exists(depth_file)) {
    depth <- read.csv(depth_file)
    colnames(depth)[1] <- "ID"
    colnames(depth)[2] <- "depth_macrofungi"

    if (is.null(depth_merged)) {
      depth_merged <- depth
    } else {
      depth_merged <- rbind(depth_merged, depth)
      rm(depth)
      gc()
    }
  } else {
    cat("No file found: ", project_name, "\n")
    not_merged_projects <- c(not_merged_projects, project_name)
  }
}

not_merged_projects
head(depth_merged)

df <- depth_merged %>% filter(depth_microfungi >= 10000)
nrow(df)
nrow(depth_merged)

write.csv(depth_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_all.csv")
write.csv(depth_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_microfungi.csv")
write.csv(depth_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_plant.csv")
write.csv(depth_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_macrofungi.csv")


# Merge raw depth with the other four depth tables, using raw as the base
depth_raw <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_raw.csv", row.names = 1)
depth_all <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_all.csv", row.names = 1)
depth_microfungi <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_microfungi.csv", row.names = 1)
depth_macrofungi <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_macrofungi.csv", row.names = 1)
depth_plant <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_plant.csv", row.names = 1)

# Select required columns from raw depth
depth_raw <- depth_raw %>% select(file, num_seqs, PRJ)
colnames(depth_raw)[1] <- "ID"
colnames(depth_raw)[2] <- "depth_raw"
depth_raw <- depth_raw %>% select(ID, PRJ, depth_raw)
head(depth_raw)

# Join depth tables
depth_merged <- depth_raw %>%
  left_join(depth_all, by = "ID") %>%
  left_join(depth_microfungi, by = "ID") %>%
  left_join(depth_macrofungi, by = "ID") %>%
  left_join(depth_plant, by = "ID") %>%
  mutate(across(everything(), ~replace_na(., 0)))

head(depth_merged)
nrow(depth_merged)

write.csv(depth_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/depth/depth_merged.csv")





# --- Merge ASV-level phyloseq objects of microfungi ---
ps_merged <- NULL  # Container for merged results
not_merged_projects <- character()

for (proj in projects) {
  project_name <- basename(proj)
  cat("Processing:", project_name, "\n")

  seqtab_file <- file.path(proj, "3_phyloseq/fungi/ASV", paste0(project_name, "_otu_fungi.csv"))
  taxa_file   <- file.path(proj, "3_phyloseq/fungi/ASV", paste0(project_name, "_taxa_fungi.csv"))
  refseq_file <- file.path(proj, "3_phyloseq/fungi/ASV", paste0(project_name, "_refseq_fungi.fasta"))

  if (file.exists(seqtab_file) && file.exists(taxa_file)) {
    seqtab <- fread(seqtab_file, header = T)
    seqtab <- as.data.frame(seqtab)
    rownames(seqtab) <- seqtab[[1]]
    seqtab <- seqtab[, -1]
    #seqtab <- t(seqtab) # Enable when merging at genus level

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
      rm(ps)         # Free memory by removing temporary objects
      gc()           # Trigger garbage collection
    }
  } else {
    cat("No file found: ", project_name, "\n")
    not_merged_projects <- c(not_merged_projects, project_name)
  }
}

not_merged_projects
ps_merged # 92050 taxa and 46598 samples

# Remove samples with 0 reads
ps_merged2 <- prune_samples(sample_sums(ps_merged) >= 1, ps_merged)
ps_merged2 <- filter_taxa(ps_merged2, function(otu) sum(otu) > 1, TRUE) # Remove ASVs that are all zeros after filtering
ps_merged2 # 92050 taxa and 46422 samples


# Write output
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV"
saveRDS(ps_merged2, file = file.path(path, "ps_merged_ASV.rds"))
write.csv(as.data.frame(otu_table(ps_merged2)), file = file.path(path, "otu_ASV.csv"))
write.csv(as.data.frame(tax_table(ps_merged2)), file = file.path(path, "tax_ASV.csv"))
writeXStringSet(ps_merged2@refseq, file = file.path(path, "refseq.fasta"), format = "fasta")


# Filter by sequencing depth to obtain ps for downstream analyses
# Remove samples with insufficient reads
ps_merged3 <- prune_samples(sample_sums(ps_merged2) >= 10000, ps_merged2)
ps_merged3 <- filter_taxa(ps_merged3, function(otu) sum(otu) > 1, TRUE) # Remove ASVs that are all zeros after filtering
ps_merged3 # 81983 taxa and 37417 samples


# Write output
path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/ASV"
saveRDS(ps_merged3, file = file.path(path, "ps_merged_ASV_filterdp.rds"))
write.csv(as.data.frame(otu_table(ps_merged3)), file = file.path(path, "otu_ASV_filterdp.csv"))
write.csv(as.data.frame(tax_table(ps_merged3)), file = file.path(path, "tax_ASV_filterdp.csv"))
writeXStringSet(ps_merged3@refseq, file = file.path(path, "refseq_filterdp.fasta"), format = "fasta")



# For the unfiltered-depth ps, build ID–PRJ mapping and compute alpha diversity
# Write the ID–PRJ mapping
id_prj <- as.matrix(sample_data(ps_merged2))
colnames(id_prj) <- c("PRJ", "ID")
head(id_prj)
nrow(id_prj)

write.csv(id_prj, "/data/wangxinyu/ITS_Public/Project/Analysis/data/meta/ID_to_PRJ.csv")

# Recompute alpha diversity consistently (Observed, Shannon, Simpson)
div <- estimate_richness(ps_merged2, split = TRUE, measures = c(c("Observed", "Shannon", "Simpson")))
head(div)
nrow(div)

write.csv(div, "/data/wangxinyu/ITS_Public/Project/Analysis/data/diversity/diversity_microfungi.csv")

# Subset diversity results to depth-filtered samples
ID_fdp <- row.names(sample_data(ps_merged3))
ID_fdp

div_fdp <- div %>%
  rownames_to_column("ID") %>%
  filter(ID %in% ID_fdp) %>%
  column_to_rownames("ID")

head(div_fdp)
nrow(div_fdp)

write.csv(div_fdp, "/data/wangxinyu/ITS_Public/Project/Analysis/data/diversity/diversity_microfungi_filterdp.csv")



# --- Merge diversity metrics of all species ---
diversity_merged <- NULL  # Container for merged results
not_merged_projects <- character()

for (proj in projects) {
  project_name <- basename(proj)
  cat("Processing:", project_name, "\n")

  diversity_file <- file.path(proj, "3_phyloseq/all", paste0(project_name, "_diversity_all.csv"))

  if (file.exists(diversity_file)) {
    diversity <- read.csv(diversity_file, row.names = 1)
    diversity$PRJ <- project_name

    if (is.null(diversity_merged)) {
      diversity_merged <- diversity
    } else {
      diversity_merged <- rbind(diversity_merged, diversity)
      rm(diversity)       # Free memory by removing temporary objects
      gc()          # Trigger garbage collection
    }
  } else {
    cat("No file found: ", project_name, "\n")
    not_merged_projects <- c(not_merged_projects, project_name)
  }
}

not_merged_projects
head(diversity_merged)

write.csv(diversity_merged, "/data/wangxinyu/ITS_Public/Project/Analysis/data/diversity/diversity_all.csv")




# --- Merge ASV-level tables from mouse and environmental datasets ---
projects <- list.dirs("/data/wangxinyu/ITS_Public/Project/Other_site/env/Done", full.names = TRUE, recursive = FALSE)
projects

ps_merged <- NULL  # Container for merged results
not_merged_projects <- character()

for (proj in projects) {
  project_name <- basename(proj)
  cat("Processing:", project_name, "\n")

  #seqtab_file <- file.path(proj, "3_phyloseq/fungi/ASV", paste0(project_name, "_otu_fungi_filterdp.csv"))
  #taxa_file   <- file.path(proj, "3_phyloseq/fungi/ASV", paste0(project_name, "_taxa_fungi_filterdp.csv"))

  seqtab_file <- file.path(proj, "3_phyloseq/fungi/genus", paste0(project_name, "_otu_fungi_gen_p0_filterdp.csv"))
  taxa_file   <- file.path(proj, "3_phyloseq/fungi/genus", paste0(project_name, "_taxa_fungi_gen_filterdp.csv"))

  if (file.exists(seqtab_file) && file.exists(taxa_file)) {
    seqtab <- fread(seqtab_file, header = T)
    seqtab <- as.data.frame(seqtab)
    rownames(seqtab) <- seqtab[[1]]
    seqtab <- seqtab[, -1]
    seqtab <- t(seqtab) # Enable when merging at genus level

    taxa <- fread(taxa_file, header = T)
    taxa <- as.data.frame(taxa)
    rownames(taxa) <- taxa[[1]]
    taxa <- taxa[, -1]

    samdf <- data.frame(cohort = project_name,
                        id = rownames(seqtab),
                        row.names = rownames(seqtab))

    ps <- phyloseq(
      otu_table(as.matrix(seqtab), taxa_are_rows = FALSE),
      sample_data(samdf),
      tax_table(as.matrix(taxa))
    )

    if (is.null(ps_merged)) {
      ps_merged <- ps
    } else {
      ps_merged <- merge_phyloseq(ps_merged, ps)
      rm(ps)         # Free memory by removing temporary objects
      gc()           # Trigger garbage collection
    }
  } else {
    cat("No file found: ", project_name, "\n")
    not_merged_projects <- c(not_merged_projects, project_name)
  }
}

ps_merged
not_merged_projects
head(otu_table(ps_merged))
head(tax_table(ps_merged))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/env"
write.csv(otu_table(ps_merged), file = file.path(path, "e_otu_ASV_filterdp.csv"))
write.csv(tax_table(ps_merged), file = file.path(path, "e_tax_ASV_filterdp.csv"))

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/env"
write.csv(otu_table(ps_merged), file = file.path(path, "e_otu_gen_p0_filterdp.csv"))
write.csv(tax_table(ps_merged), file = file.path(path, "e_tax_gen_p0_filterdp.csv"))



# --- Merge datasets from other sites and human ---
# human
otu_human <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/genus/otu_gen_p0_filterdp.csv")
otu_human <- t(otu_human)
otu_human[1:3, 1:3]

tax_human <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/genus/tax_gen_p0_filterdp.csv", row.names = 1)
tax_human[1:3, 1:3]

meta_human <- data.frame(site = "human_gut",
                         id = rownames(otu_human),
                         row.names = rownames(otu_human))
head(meta_human)

ps_human <- phyloseq(otu_table(as.matrix(otu_human), taxa_are_rows = FALSE),
                     tax_table(as.matrix(tax_human)),
                     sample_data(meta_human))

ps_human

# mouse
otu_mouse <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/mouse/m_otu_gen_p0_filterdp.csv")
otu_mouse[1:3, 1:3]

tax_mouse <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/mouse/m_tax_gen_p0_filterdp.csv", row.names = 1)
tax_mouse[1:3, 1:3]

meta_mouse <- data.frame(site = "mouse_gut",
                   id = rownames(otu_mouse),
                   row.names = rownames(otu_mouse))
head(meta_mouse)

ps_mouse <- phyloseq(otu_table(as.matrix(otu_mouse), taxa_are_rows = FALSE),
                     tax_table(as.matrix(tax_mouse)),
                     sample_data(meta_mouse))

ps_mouse

# Environment
otu_env <- read_otu("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/env/e_otu_gen_p0_filterdp.csv")
otu_env[1:3, 1:3]

tax_env <- read.csv("/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/env/e_tax_gen_p0_filterdp.csv", row.names = 1)
tax_env[1:3, 1:3]

meta_env <- data.frame(site = "environment",
                   id = rownames(otu_env),
                   row.names = rownames(otu_env))
head(meta_env)

ps_env <- phyloseq(otu_table(as.matrix(otu_env), taxa_are_rows = FALSE),
                     tax_table(as.matrix(tax_env)),
                     sample_data(meta_env))

ps_env


# merge
ps_merged <- merge_phyloseq(ps_human, ps_mouse, ps_env)
ps_merged

meta_merged <- as.matrix(sample_data(ps_merged))
head(meta_merged)
class(meta_merged)

path <- "/data/wangxinyu/ITS_Public/Project/Analysis/data/other_site/merged"
write.csv(otu_table(ps_merged), file = file.path(path, "otu_gen_p0_filterdp_merged.csv"))
write.csv(tax_table(ps_merged), file = file.path(path, "tax_gen_p0_filterdp_merged.csv"))
write.csv(meta_merged, file = file.path(path, "meta_gen_p0_filterdp_merged.csv"))
