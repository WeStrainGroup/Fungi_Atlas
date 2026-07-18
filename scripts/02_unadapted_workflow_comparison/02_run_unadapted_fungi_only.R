#!/usr/bin/env Rscript
# Purpose: Run the primary fixed-trimming DADA2 workflow with a fungi-only UNITE database.

suppressMessages({
  library(R.utils)
  library(tidyverse)
  library(data.table)
  library(fs)
  library(dada2)
  library(Biostrings)
  library(phyloseq)
  library(microbiome)
  library(vegan)
})

project0 <- "WeP1"
path_in <- "/data/wangxinyu/Fungi_Atlas/Data/WeGut_raw/WeP1/dada/raw"
path_out <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted1/WeP1"
pattern_F <- ".1.fq.gz"
pattern_R <- ".2.fq.gz"
trimLeft_F <- 26 #HMP 55, CHGM 28, WeP1 26
trimLeft_R <- 26 #HMP 30, CHGM 26, WeP1 26

# project0 <- "HMP"
# path_in <- "/data/wangxinyu/Fungi_Atlas/Data/Public_raw/PE/PRJNA356769b_ITS2"
# path_out <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted1/HMP"
# pattern_F <- ".1.fastq.gz"
# pattern_R <- ".2.fastq.gz"
# trimLeft_F <- 55 #HMP 55, CHGM 28, WeP1 26
# trimLeft_R <- 30 #HMP 30, CHGM 26, WeP1 26

# project0 <- "CHGM"
# path_in <- "/data/wangxinyu/Fungi_Atlas/Data/Public_raw/PE/PRJCA010668_ITS1"
# path_out <- "/data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/unadapted1/CHGM"
# pattern_F <- ".1.fastq.gz"
# pattern_R <- ".2.fastq.gz"
# trimLeft_F <- 28 #HMP 55, CHGM 28, WeP1 26
# trimLeft_R <- 26 #HMP 30, CHGM 26, WeP1 26



project <- paste0(project0, "_")
if(!dir.exists(path_out)) {dir.create(path_out)}
seq_type <- "PE"
maxEE_F <- 5
maxEE_R <- 5

# use fungi only database
ref <- "/data/wangxinyu/Ref/ITS_anno/Unite/sh_general_release_19.02.2025/sh_general_release_dynamic_19.02.2025.fasta"
minBoot <- 50
filter_abundance <- 1/10000
filter_depth <- 10000
thread <- 16


message("\n=================================================")
message("\nRunning with the following parameters — make sure they're correct!\n")
message(paste0("Project name: ", project0))
message(paste0("Input directory: ", path_in))
message(paste0("Output directory: ", path_out))
message(paste0("Sequencing type: ", seq_type))
message(paste0("Forward reads pattern: ", pattern_F))
message(paste0("Reverse reads pattern: ", pattern_R))
message(paste0("Maximum expected errors (maxEE) for forward reads filtering: ", maxEE_F))
message(paste0("Maximum expected errors (maxEE) for reverse reads filtering: ", maxEE_R))
message(paste0("Forward reads: trimming ", trimLeft_F, " nucleotide(s) from the start (trimLeft)."))
message(paste0("Reverse reads: trimming ", trimLeft_R, " nucleotide(s) from the start (trimLeft)."))
message(paste0("Reference database: ", basename(ref)))
message(paste0("Minimum bootstrap value: ", minBoot))
message(paste0("Relative abundance filter threshold: ", filter_abundance))
message(paste0("Sequencing depth filter threshold: ", filter_depth))
message(paste0("Number of threads: ", thread))
message("\n=================================================")


# --- 0.2 R package dependency validation ---
required_packages <- c(
  "dada2", "Biostrings", "tidyverse", "data.table",
  "phyloseq", "microbiome", "vegan", "fs", "R.utils"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "not found. Install via Conda or R."))
  }
}

# --- 0.3 External tools configuration and validation ---
# Get tool paths dynamically (via Conda environment PATH)
seqkit <- Sys.which("seqkit")
itsx <- Sys.which("ITSx")
vsearch <- Sys.which("vsearch")

# Set ITSx environment variable
itsx_bin <- dirname(itsx)
Sys.setenv(PATH = paste(Sys.getenv("PATH"), itsx_bin, sep = ":"))

# Check if required tools exist
check_tool <- function(tool, name) {
  if (tool == "") {
    stop(paste0("Error: ", name, " not found in PATH. Ensure it is installed via Conda."))
  }
}
check_tool(seqkit, "seqkit")
check_tool(itsx, "ITSx")
check_tool(vsearch, "vsearch")

# --- 0.4 Function definitions for quality control and visualization ---
# Define function for FASTQ quality assessment and statistics
check_fq_stats <- function(F, R, batch_size = 1000) {
    title <- paste0(project, "fq_stats_", sub(".*_", "", substitute(F)))
    p_in <- path_common(path_dir(c(F, R)))  # Find common path
    all_files <- list.files(p_in, full.names = TRUE)
    p_out <- file.path(path_check, paste0(title, ".tsv"))

    # Batch division
    file_batches <- split(all_files, ceiling(seq_along(all_files) / batch_size))

    # Temporary result file list
    temp_outputs <- character(length(file_batches))

    for (i in seq_along(file_batches)) {
        batch_files <- file_batches[[i]]
        temp_out <- tempfile(pattern = paste0("fq_stats_batch_", i), fileext = ".tsv")
        temp_outputs[i] <- temp_out

        arg_seqkit <- c("stats", "-a", "-b", "-e", "-T", "-j", thread, batch_files)
        out <- system2(seqkit, args = arg_seqkit, stdout = TRUE)
        writeLines(out, temp_out)
    }

    # Merge all results and save
    df_list <- lapply(temp_outputs, read.delim, header = TRUE)
    fa_raw_check <- do.call(rbind, df_list)
    write.table(fa_raw_check, file = p_out, sep = "\t", row.names = FALSE, quote = FALSE)

    # Plot
    his1 <- ggplot(fa_raw_check, aes(x = num_seqs)) +
        geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))),
                       fill = "grey", color = "black", alpha = 0.75) +
        labs(title = title,
             x = "Num_seqs of each sample",
             y = "Proportion") +
        theme_bw()

    ggsave(file.path(path_check, paste0(title, "1.jpg")), his1, width = 6, height = 6)

    his2 <- ggplot(fa_raw_check, aes(x = avg_len)) +
            geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))), fill = "grey", color = "black", alpha = 0.75) +
            labs(title = title,
                x = "Avg_len of each sample",
                y = "Proportion") +
            theme_bw()
      ggsave(file.path(path_check,  paste0(title, "2.jpg")), his2, width = 6, height = 6)

    his3 <- ggplot(fa_raw_check, aes(x = AvgQual)) +
            geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))), fill = "grey", color = "black", alpha = 0.75) +
            labs(title = title,
                x = "Avg_qual of each sample",
                y = "Proportion") +
            theme_bw()
    ggsave(file.path(path_check, paste0(title, "3.jpg")), his3, width = 6, height = 6)
}

# Define function for quality profile visualization and export
if (seq_type == 'PE') {
    save_quality_plot <- function(fq_F, fq_R) {
    title_F <- paste0(project, "quality_", substitute(fq_F))
    title_R <- paste0(project, "quality_", substitute(fq_R))

    p_quality_R1 <- plotQualityProfile(fq_F[1:9]) + labs(title = title_F)
    p_quality_R2 <- plotQualityProfile(fq_R[1:9]) + labs(title = title_R)

    ggsave(file.path(path_check, paste0(title_F, ".jpg")), p_quality_R1, width = 6, height = 6)
    ggsave(file.path(path_check, paste0(title_R, ".jpg")), p_quality_R2, width = 6, height = 6)
  }
} else {
    save_quality_plot <- function(fq) {
    title <- paste0(project, "quality_", substitute(fq))
    p_quality <- plotQualityProfile(fq[1:9]) + labs(title = title)
    ggsave(file.path(path_check, paste0(title, ".jpg")), p_quality, width = 6, height = 6)
  }
}




# Echo both command and output
options(echo = TRUE)

# --- 1 DADA2-based sequence processing pipeline ---
path1 <- file.path(path_out, "1_dada2")
if(!dir.exists(path1)) {dir.create(path1)}

  # Input FASTQ files
  fqF_raw <- sort(list.files(path_in, pattern = pattern_F, full.names = TRUE))
  fqR_raw <- sort(list.files(path_in, pattern = pattern_R, full.names = TRUE))

  length(fqF_raw)
  head(fqF_raw)
  length(fqR_raw)
  head(fqR_raw)

  path_check <- file.path(path1, "check")
  if(!dir.exists(path_check)) {dir.create(path_check)}

  check_fq_stats(fqF_raw, fqR_raw)
  save_quality_plot(fqF_raw, fqR_raw)

  # --- 1.1 Quality filtering and trimming using DADA2 ---
  path_filtered <- file.path(path1, "filtered")
  if(!dir.exists(path_filtered)) dir.create(path_filtered)

  fqF_filtered <- file.path(path_filtered, basename(fqF_raw))
  fqR_filtered <- file.path(path_filtered, basename(fqR_raw))

  out <- filterAndTrim(fqF_raw, fqF_filtered, fqR_raw, fqR_filtered,
                      trimLeft = c(trimLeft_F, trimLeft_R),
                      truncQ = 2, maxN = 0, maxEE = c(maxEE_F, maxEE_R),
                      rm.phix = TRUE, compress = T, multithread = thread, verbose = TRUE)

  rownames(out) <- gsub(pattern_F, "", rownames(out))
  head(out)

  # Keep only FASTQ files with reads
  fqF_filtered <- sort(list.files(path_filtered, pattern = pattern_F, full.names = TRUE))
  fqR_filtered <- sort(list.files(path_filtered, pattern = pattern_R, full.names = TRUE))

  # Check
  check_fq_stats(fqF_filtered, fqR_filtered)
  save_quality_plot(fqF_filtered, fqR_filtered)

  # --- 1.2 DADA2 error rate estimation and denoising ---
  # Learn the Error Rates
  errF <- learnErrors(fqF_filtered, nbases = 1e9, randomize = TRUE, multithread = thread)
  errR <- learnErrors(fqR_filtered, nbases = 1e9, randomize = TRUE, multithread = thread)

  # Plot error rates
  p_error_F <- plotErrors(errF, nominalQ = TRUE) + labs(title = paste0(project, "error_F"))
  p_error_R <- plotErrors(errR, nominalQ = TRUE) + labs(title = paste0(project, "error_R"))

  ggsave(file.path(path_check, paste0(project, "error_F.png")), p_error_F, width = 6, height = 6)
  ggsave(file.path(path_check, paste0(project, "error_R.png")), p_error_R, width = 6, height = 6)

  # Sample inference to the de-replicated data
  dadaFs <- dada(fqF_filtered, err = errF, multithread = thread)
  dadaRs <- dada(fqR_filtered, err = errR, multithread = thread)

  dadaFs[[1]]
  dadaRs[[1]]

  # --- 1.3 Paired-end read merging and ASV table construction ---
  # Merge paired reads
  mergers <- mergePairs(dadaFs, fqF_filtered, dadaRs, fqR_filtered, verbose = TRUE)
  head(mergers[[1]])

  # Construct sequence (ASV) table
  seqtab <- makeSequenceTable(mergers)
  rownames(seqtab) <- gsub(pattern_F, "", rownames(seqtab))

  dim(seqtab)
  table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

  # Write output
  path_result <- file.path(path1, "result")
  if(!dir.exists(path_result)) dir.create(path_result)

  write.csv(seqtab, file = file.path(path_result, paste0(project, "seqtab_raw.csv")))



# --- 1.5 Chimera removal and final ASV table generation ---
# Use a larger minFoldParentOverabundance value for PB data
if (seq_type == 'PB') {
  minFoldParentOverAbundance = 3.5
} else {
  minFoldParentOverAbundance = 2
}

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", minFoldParentOverAbundance = minFoldParentOverAbundance,
       multithread = thread, verbose = TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
table(nchar(getSequences(seqtab.nochim)))

# Update sample names
row.names(seqtab.nochim) <- gsub(pattern_F, "", row.names(seqtab.nochim))
head(row.names(seqtab.nochim))

# Write output
write.csv(seqtab.nochim, file = file.path(path_result, paste0(project, "seqtab.csv")))


# --- 1.6 Taxonomic assignment ---
seqs <- colnames(seqtab.nochim)
length(seqs)
head(seqs)

batch_size <- thread * 100
batch_size
num_batches <- ceiling(length(seqs) / batch_size)
num_batches

taxa_list <- list()
boot_list <- list()

for (i in 1:num_batches) {
  up <- min(i * batch_size, length(seqs))
  down <- (i - 1) * batch_size + 1
  seqs_batch <- seqs[down:up] # Select sequences for this batch
  taxa_batch <- assignTaxonomy(seqs_batch, ref, minBoot = minBoot, outputBootstraps = TRUE, multithread = thread, tryRC = TRUE, verbose = TRUE)

  taxa_list[[i]] <- taxa_batch[["tax"]]
  boot_list[[i]] <- taxa_batch[["boot"]]

  rm(seqs_batch, taxa_batch) # Free up memory
  cat(paste0("\n", "-- batch", i, " / ", num_batches,": seq", down, "-", up ," was down! --", "\n"))
}

taxa <- do.call(rbind, taxa_list)
boot <- do.call(rbind, boot_list)

nrow(taxa)
nrow(boot)
sum(is.na(taxa[, "Species"]))
sum(is.na(taxa[, "Species"]))/nrow(taxa)
sum(is.na(taxa[, "Genus"]))
sum(is.na(taxa[, "Genus"]))/nrow(taxa)

# Write output
write.csv(taxa, file = file.path(path_result, paste0(project, "taxa.csv")))
write.csv(boot, file = file.path(path_result, paste0(project, "taxa_boot.csv")))



# --- 3 Phyloseq-based community analysis pipeline ---
path3 <- file.path(path_out, "3_phyloseq")
if(!dir.exists(path3)) {dir.create(path3)}

# Construct the phyloseq subject
otu_table <- phyloseq::otu_table(as.matrix(seqtab.nochim), taxa_are_rows = FALSE)
otu_table[1:3, 1:3]

sam_data <- data.frame(row.names = rownames(otu_table), project = rep(project0, nrow(otu_table)))
sam_data <- sample_data(sam_data)

#taxa <- read.csv("/data/wangxinyu/Q_ITS_Error/Anno_classic/PRJNA356769b_ITS2_classic/1_dada2/result/PRJNA356769b_ITS2_classic_taxa.csv", row.names = 1)
tax_table <- phyloseq::tax_table(as.matrix(taxa))
tax_table[1:3, 1:3]


ps <- phyloseq(otu_table,
               sam_data,
               tax_table)
ps
microbiome::summarize_phyloseq(ps)

# --- 3.1 Global community analysis and basic filtering ---
path_all <- file.path(path3, "all")
if(!dir.exists(path_all)) dir.create(path_all)

# --- 3.1.1 Do general filteration  ---
# Remove ASVs with unclear phylum annotation
ps_filter <- subset_taxa(ps, Kingdom != "k__Eukaryota_kgd_Incertae_sedis")

ntaxa(ps)
ntaxa(ps_filter)

# Set ASV counts to 0 in samples where relative abundance is below filter_abundance
otu <- ps_filter@otu_table
otu[otu / rowSums(otu) < filter_abundance] <- 0
otu_table(ps_filter) <- otu # Return otu table to ps
ps_filter <- filter_taxa(ps_filter, function(otu) sum(otu) > 1, TRUE) # Remove ASVs that are all zeros after filtering
ps_filter <- prune_samples(sample_sums(ps_filter) > 0, ps_filter) # Remove samples without any reads

ntaxa(ps)
ntaxa(ps_filter)
nsamples(ps)
nsamples(ps_filter)

# --- 3.1.2 Write output  ---
path_all
saveRDS(ps_filter, file.path(path_all, paste0(project, "ps_all.rds")))
write.csv(readcount(ps_filter), file.path(path_all, paste0(project, "readcount_all.csv")))
write.csv(ps_filter@tax_table, file.path(path_all, paste0(project, "taxa_all.csv")))
write.csv(ps_filter@otu_table, file.path(path_all, paste0(project,"otu_all.csv")))

# Alpha diversity
try({diversity_all <- estimate_richness(ps_filter, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))
     write.csv(diversity_all, file.path(path_all, paste0(project, "diversity_all.csv")))
     }, silent = F)


options(echo = FALSE)

message("===============================================\n")
message("All done!\n")
message(paste0("End at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n"))
message("===============================================\n")
