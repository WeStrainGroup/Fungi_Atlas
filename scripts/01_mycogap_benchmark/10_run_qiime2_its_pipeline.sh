#!/bin/bash
# Purpose: Process mock reads with QIIME2, Cutadapt, ITSxpress, DADA2, VSEARCH clustering, and UNITE classification.


prj="mock"
threads=16
path_in="/home/wangxy/ITS_benchmark/data/fq"
path_out="/home/wangxy/ITS_benchmark/process/qiime2"
path_out1="${path_out}/process"
path_out2="${path_out}/result"

mkdir -p $path_out1 $path_out2


echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${path_out1}/${prj}_list.tsv"

for f in "$path_in"/*_R1.fastq.gz; do
    s=$(basename "$f" _R1.fastq.gz)
    echo -e "$s\t$(realpath "$f")\t$(realpath "$path_in/${s}_R2.fastq.gz")"
done >> "${path_out1}/${prj}_list.tsv"


qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path "${path_out1}/${prj}_list.tsv" \
--output-path "${path_out1}/${prj}_raw.qza" \
--input-format PairedEndFastqManifestPhred33V2


qiime cutadapt trim-paired \
--i-demultiplexed-sequences "${path_out1}/${prj}_raw.qza" \
--o-trimmed-sequences "${path_out1}/${prj}_cut.qza" \
--p-front-f GCATCGATGAAGAACGCAGC \
--p-front-r TCCTCCGCTTATTGATATGC \
--p-discard-untrimmed \
--p-cores $threads


qiime itsxpress trim-pair-output-unmerged \
--i-per-sample-sequences "${path_out1}/${prj}_cut.qza" \
--o-trimmed "${path_out1}/${prj}_its.qza" \
--p-region ITS2 \
--p-taxa ALL \
--p-cluster-id 1.0 \
--p-threads $threads


qiime dada2 denoise-paired \
--i-demultiplexed-seqs "${path_out1}/${prj}_its.qza" \
--p-trunc-len-f 0 \
--p-trunc-len-r 0 \
--p-max-ee-f 5 \
--p-max-ee-r 5 \
--p-trim-overhang TRUE \
--p-n-reads-learn 1000000000 \
--p-n-threads $threads \
--o-representative-sequences "${path_out1}/${prj}_seq.qza" \
--o-table "${path_out1}/${prj}_seqtab.qza" \
--o-denoising-stats "${path_out1}/${prj}_denoising_stats.qza" \
--o-base-transition-stats "${path_out1}/${prj}_transition_stats.qza"

# post-cluster
qiime vsearch cluster-features-de-novo \
--i-sequences "${path_out1}/${prj}_seq.qza" \
--i-table "${path_out1}/${prj}_seqtab.qza" \
--o-clustered-sequences "${path_out1}/${prj}_seq_cluster.qza" \
--o-clustered-table "${path_out1}/${prj}_seqtab_cluster.qza" \
--p-perc-identity 0.985 \
--p-strand both \
--p-threads $threads


qiime feature-classifier classify-sklearn \
--i-reads "${path_out1}/${prj}_seq_cluster.qza" \
--i-classifier "/home/wangxy/ITS_benchmark/ref/qiime2/v2026.1/unite_ver2025-02-19_dynamic_eukaryotes-Q2-2026.1.qza" \
--p-n-jobs $threads \
--o-classification "${path_out1}/${prj}_taxtab_cluster.qza"


qiime tools export \
--input-path "${path_out1}/${prj}_seqtab_cluster.qza" \
--output-path "${path_out2}/"

qiime tools export \
--input-path "${path_out1}/${prj}_seq_cluster.qza" \
--output-path "${path_out2}/"

qiime tools export \
--input-path "${path_out1}/${prj}_taxtab_cluster.qza" \
--output-path "${path_out2}/"
