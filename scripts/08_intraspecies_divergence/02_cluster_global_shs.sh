#!/bin/bash
# Purpose: Cluster pooled ITS1 and ITS2 ASVs separately at 98.5% identity to define global SHs.


vsearch \
--cluster_fast /data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_refseq_fdp_ITS1.fasta \
--centroids /data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS1_v98.5.fasta \
--uc /data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS1_v98.5.uc \
--id 0.985 \
--strand both \
--threads 96


vsearch \
--cluster_fast /data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/raw/unc_refseq_fdp_ITS2.fasta \
--centroids /data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS2_v98.5.fasta \
--uc /data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance/vsearch/unc_refseq_fdp_ITS2_v98.5.uc \
--id 0.985 \
--strand both \
--threads 96
