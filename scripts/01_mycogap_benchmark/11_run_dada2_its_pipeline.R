#!/usr/bin/env Rscript
# Purpose: Process mock reads with the DADA2_ITS comparison workflow and save a phyloseq object.



library(dada2)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(tidyverse)
library(microbiome)

options(echo = TRUE)

thread <- 16

# input
path_in <- "/home/wangxy/ITS_benchmark/data/fq"
path_out <- "/home/wangxy/ITS_benchmark/process/dada2"
list.files(path_in, pattern = ".fastq.gz")


fnFs <- sort(list.files(path_in, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path_in, pattern = "_R2.fastq.gz", full.names = TRUE))
fnFs
fnRs


# identify primers
FWD <- "GCATCGATGAAGAACGCAGC"
REV <- "TCCTCCGCTTATTGATATGC"

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# pre filter and trim
fnFs.filtN <- file.path(path_out, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path_out, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = thread)


# count primer
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


# remove primer
cutadapt <- "/home/wangxy/miniconda3/envs/mycogap/bin/cutadapt"
system2(cutadapt, args = "--version")


path.cut <- file.path(path_out, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "--cores", thread,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))


# main filter and trim
filtFs <- file.path(path_out, "filtered", basename(cutFs))
filtRs <- file.path(path_out, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN = 0, maxEE = c(5, 5), truncQ = 2, minLen = 50, rm.phix = TRUE,
                     compress = TRUE, multithread = thread)
head(out)

# check
p_quality_R1 <- plotQualityProfile(filtFs[1:9]) + labs(title = "quality_F")
p_quality_R2 <- plotQualityProfile(filtRs[1:9]) + labs(title = "quality_R")

ggsave(file.path(path_out, "quality_F.jpg"), p_quality_R1, width = 6, height = 6)
ggsave(file.path(path_out, "quality_R.jpg"), p_quality_R2, width = 6, height = 6)

# lean the error rate
errF <- learnErrors(filtFs, nbases = 1e9, randomize = TRUE, multithread = thread)
errR <- learnErrors(filtRs, nbases = 1e9, randomize = TRUE, multithread = thread)

# check
p_error_F <- plotErrors(errF, nominalQ = TRUE) + labs(title = "error_F")
p_error_R <- plotErrors(errR, nominalQ = TRUE) + labs(title = "error_R")

ggsave(file.path(path_out, "error_F.png"), p_error_F, width = 6, height = 6)
ggsave(file.path(path_out, "error_R.png"), p_error_R, width = 6, height = 6)

# sample interence
dadaFs <- dada(filtFs, err = errF, multithread = thread)
dadaRs <- dada(filtRs, err = errR, multithread = thread)

# merge pair
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

# construct feature table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = thread, verbose = TRUE)
rownames(seqtab.nochim) <- gsub("_R1.fastq.gz", "", rownames(seqtab.nochim))

seqtab.nochim[1:3, 1:3]
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))


# track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- gsub("_R1.fastq.gz", "", rownames(track))

head(track)

# assign taxonomy
unite.ref <- "/home/wangxy/ITS_benchmark/ref/blast/sh_general_release_dynamic_all_19.02.2025.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = thread, tryRC = TRUE)
head(taxa)

# construct the phyloseq subject
otu_table <- phyloseq::otu_table(seqtab.nochim, taxa_are_rows = FALSE)
sam_data <- sample_data(data.frame(row.names = rownames(otu_table), project = rep("mock", nrow(otu_table))))
tax_table <- phyloseq::tax_table(as.matrix(taxa))

ps <- phyloseq(otu_table,
               sam_data,
               tax_table)
ps

# add ASV tag
ps <- add_refseq(ps, tag = "ASV")
taxa_names(ps)
ps

# write out
saveRDS(ps, file.path(path_out, "mock_dada2_its.rds"))
