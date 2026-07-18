# Purpose: Curate 100%-matched references and generate mutation-aware mock sequences and log-normal abundance profiles.


library(tidyverse)
library(Biostrings)



blast_res <- read.table("/home/wangxy/ITS_benchmark/data/blast/refseq_fun_obipcr_its2_v97_blast.tsv", header = 1)
head(blast_res)


blast_res_filter <- blast_res %>% filter(pident == 100)

nrow(blast_res)
nrow(blast_res_filter)

write.csv(blast_res_filter, "/home/wangxy/ITS_benchmark/data/blast/refseq_fun_obipcr_its2_v97_blast_filter.csv")


fa <- readDNAStringSet("/home/wangxy/ITS_benchmark/data/refseq_fun_obipcr_its2_v97.fa", format = "fasta")
fa_filter <- fa[names(fa) %in% blast_res_filter$qseqid]

fa
fa_filter

writeXStringSet(fa_filter, "/home/wangxy/ITS_benchmark/data/refseq_fun_obipcr_its2_v97_filter.fa")



# Run 06_sample_mock_communities.sh at this checkpoint, then continue below.




df <- read.csv("/home/wangxy/ITS_benchmark/data/data_s1_iscience_2023.csv")
head(df)


df_filter <- df %>%
    select(Assembly.Reference, ITS.Copies, Pairwise..Identity) %>%
    filter(ITS.Copies >= 3) %>%
    filter(Pairwise..Identity >= 97) %>%
    arrange(Pairwise..Identity) %>%
    mutate(group = case_when(
      Pairwise..Identity == 100 ~ 100,
      Pairwise..Identity >= 99 & Pairwise..Identity < 100 ~ 99,
      Pairwise..Identity >= 98 & Pairwise..Identity < 99 ~ 98,
      Pairwise..Identity >= 97 & Pairwise..Identity < 98 ~ 97
    )
  )

nrow(df) # 2414
nrow(df_filter) # 412
head(df_filter)
table(df_filter$group)

dis <- round(prop.table(table(df_filter$group)) * 100, 0)
dis <- round(prop.table(table(df_filter$group)) * 100, 2)
dis


round(weighted.mean(as.numeric(names(dis)), as.vector(dis)), 1) #  99





# conda activate genometools
indir <- "/home/wangxy/ITS_benchmark/data/sampling"
outdir <- "/home/wangxy/ITS_benchmark/data/fq"
dir.create(outdir, showWarnings = FALSE)

fa_inputs <- list.files(indir, pattern = "\\.fa$", full.names = TRUE)
fa_inputs

for(infile in fa_inputs){

  cat("doing:", infile, "\n")

  prefix <- tools::file_path_sans_ext(basename(infile))

  mutdir <- file.path("/home/wangxy/ITS_benchmark/data/seqmutate", prefix)
  dir.create(mutdir, showWarnings = FALSE, recursive = TRUE)

  fa <- readDNAStringSet(infile)

  mut_times <- rep(1, length(fa))

  mut_rates <- sample(
    c(3, 2, 1, 0),
    size = length(fa),
    replace = TRUE,
    prob = c(0.05, 0.13, 0.55, 0.27)
  )


  fa_files <- c()


  for(i in seq_along(fa)){

    seq_name <- paste0(prefix, "_", names(fa)[i])

    seq_list <- list()
    tmp_files <- c()


    original_file <- file.path(mutdir, paste0(seq_name, "_mut0.fa"))
    writeXStringSet(fa[i], original_file)

    seq0 <- readDNAStringSet(original_file)
    names(seq0) <- paste0(seq_name, "_mut0_r0")

    seq_list[[1]] <- seq0
    tmp_files <- c(tmp_files, original_file)


    for(j in seq_len(mut_times[i])){

      out_file <- file.path(mutdir, paste0(seq_name, "_mut", j, ".fa"))

      rate_j <- mut_rates[i]

      cmd <- sprintf(
        'gt seqmutate -force -rate %d -o "%s" "%s"',
        rate_j, out_file, original_file
      )
      system(cmd)

      mut_seq <- readDNAStringSet(out_file)
      names(mut_seq) <- paste0(seq_name, "_mut", j, "_r", rate_j)

      seq_list[[j + 1]] <- mut_seq
      tmp_files <- c(tmp_files, out_file)
    }


    all_seq <- do.call(c, seq_list)

    final_file <- file.path(mutdir, paste0(seq_name, "_mut.fa"))
    writeXStringSet(all_seq, final_file)

    fa_files <- c(fa_files, final_file)


    file.remove(tmp_files)
    file.remove(list.files(mutdir, pattern = "_mut0", full.names = TRUE))

    cat("done:", seq_name,
        "mut_time =", mut_times[i],
        "mut_rate =", mut_rates[i], "\n")
  }


  all_fa <- do.call(c, lapply(fa_files, readDNAStringSet))
  out_fa <- file.path(outdir, paste0(prefix, "_sequence.fa"))
  writeXStringSet(all_fa, out_fa)




  sdlog <- runif(1, 0.5, 1.5)
  taxa_abundance <- rlnorm(length(fa_files), meanlog = 2, sdlog = sdlog)
  taxa_abundance <- taxa_abundance / sum(taxa_abundance)

  abundance_list <- list()

  for(i in seq_along(fa_files)){

    seqs <- readDNAStringSet(fa_files[i])
    n_copy <- length(seqs)

    per_copy <- taxa_abundance[i] / n_copy

    df <- data.frame(
      sequence = names(seqs),
      abundance = rep(per_copy, n_copy),
      stringsAsFactors = FALSE
    )

    abundance_list[[i]] <- df
  }

  abundance_table <- do.call(rbind, abundance_list)

  out_abun <- file.path(outdir, paste0(prefix, "_abundance.txt"))

  write.table(
    abundance_table,
    file = out_abun,
    sep = " ",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  cat("done:", prefix, "\n")
}







# others



# rtrunc <- function(n, mean=100, sd=60, min=10, max=250){
#  x <- rnorm(n, mean, sd)
#  x[x < min] <- min
#  x[x > max] <- max
#  round(x)
#}

# mut_times <- rtrunc(length(fa))
# mut_times



# taxa_abundance <- rep(1, length(fa_files))
