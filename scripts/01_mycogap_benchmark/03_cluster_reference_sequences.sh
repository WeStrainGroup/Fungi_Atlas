#!/bin/bash
# Purpose: Cluster validated ITS2 references at 97% identity and retain sequences 325-575 bp long.



vsearch \
--cluster_fast /home/wangxy/ITS_benchmark/refseq_ITS/refseq_fun_obipcr_its2.fa \
--centroids /home/wangxy/ITS_benchmark/refseq_ITS/refseq_fun_obipcr_its2_v97.fa \
--id 0.97 \
--minseqlength 325 \
--maxseqlength 575
