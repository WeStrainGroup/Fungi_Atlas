#!/bin/bash
# Purpose: Dereplicate mutation-aware mock truth, build sample-specific databases, and BLAST recovered MycoGAP features.


for f in /home/wangxy/ITS_benchmark/data/fq/*.fa; do
  base=$(basename "$f" .fa)

  echo "Processing $base ..."

  # 1. derep
  vsearch \
    --derep_fulllength "$f" \
    --output "/home/wangxy/ITS_fp_check/data/feature_derep/${base}_derep.fa" \
    --sizeout \
    --uc "/home/wangxy/ITS_fp_check/data/feature_derep/${base}_derep.uc"


  dbdir="/home/wangxy/ITS_fp_check/data/blast_db/${base}_derep_blast_db"
  mkdir -p "$dbdir"

  makeblastdb \
    -in "/home/wangxy/ITS_fp_check/data/feature_derep/${base}_derep.fa" \
    -dbtype nucl \
    -out "${dbdir}/${base}_derep_blast_db"


  outdir="/home/wangxy/ITS_fp_check/data/blast_res"
  mkdir -p "$outdir"

  blastn \
    -query "/home/wangxy/ITS_fp_check/data/feature_mycogap/${base}.fa" \
    -db "${dbdir}/${base}_derep_blast_db" \
    -out "${outdir}/${base}_blast0.tsv" \
    -outfmt '6 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore' \
    -max_target_seqs 1 \
    -num_threads 16


  (
    echo -e "qseqid\tsseqid\tpident\tlength\tqcovs\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue\tbitscore"
    cat "${outdir}/${base}_blast0.tsv"
  ) > "${outdir}/${base}_blast.tsv"


  rm "${outdir}/${base}_blast0.tsv"

done
