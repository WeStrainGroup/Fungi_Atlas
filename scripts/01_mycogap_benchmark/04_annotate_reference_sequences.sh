#!/bin/bash
# Purpose: Annotate clustered ITS2 references against UNITE with BLASTN for ground-truth curation.



#makeblastdb \
#-in /home/wangxy/ITS_benchmark/ref/blast/sh_general_release_dynamic_all_19.02.2025.fasta \
#-dbtype nucl \
#-out /home/wangxy/ITS_benchmark/ref/blast/sh_general_release_dynamic_all_19.02.2025_db

# run blast
blastn \
-query /home/wangxy/ITS_benchmark/data/refseq_fun_obipcr_its2_v97.fa \
-db /home/wangxy/ITS_benchmark/ref/blast/sh_general_release_dynamic_all_19.02.2025_db \
-out /home/wangxy/ITS_benchmark/data/blast/refseq_fun_obipcr_its2_v97_blast0.tsv \
-outfmt '6 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore' \
-max_target_seqs 1 \
-num_threads 64


(echo -e "qseqid\tsseqid\tpident\tlength\tqcovs\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue\tbitscore"; cat /home/wangxy/ITS_benchmark/data/blast/refseq_fun_obipcr_its2_v97_blast0.tsv) > /home/wangxy/ITS_benchmark/data/blast/refseq_fun_obipcr_its2_v97_blast.tsv
rm /home/wangxy/ITS_benchmark/data/blast/refseq_fun_obipcr_its2_v97_blast0.tsv
