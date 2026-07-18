#!/bin/bash
# Purpose: BLAST CGF genome-derived ITS1 and ITS2 sequences against the non-redundant atlas ITS database.




makeblastdb \
-in /data/wangxinyu/Fungi_Atlas/Analysis/Main/data/ASV/refseq_filterdp_derep.fasta \
-dbtype nucl \
-out /data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/fungi_atlas_blast_db/blast/fungi_atlas_blast_db


# run blast
cd /data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/blast

# ITS1
blastn \
-query /data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/ITS_CGF/CGF_ITS1.fasta \
-db fungi_atlas_blast_db/fungi_atlas_blast_db \
-out CGF_ITS1_blast0.tsv \
-outfmt '6 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore' \
-max_target_seqs 1 \
-num_threads 16

(echo -e "qseqid\tsseqid\tpident\tlength\tqcovs\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue\tbitscore"; cat CGF_ITS1_blast0.tsv) > CGF_ITS1_blast.tsv
rm  CGF_ITS1_blast0.tsv

# ITS2
blastn \
-query /data/wangxinyu/Fungi_Atlas/Analysis/Revision/compare_to_CGF/ITS_CGF/CGF_ITS2.fasta \
-db fungi_atlas_blast_db/fungi_atlas_blast_db \
-out CGF_ITS2_blast0.tsv \
-outfmt '6 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore' \
-max_target_seqs 1 \
-num_threads 16

(echo -e "qseqid\tsseqid\tpident\tlength\tqcovs\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue\tbitscore"; cat CGF_ITS2_blast0.tsv) > CGF_ITS2_blast.tsv
rm  CGF_ITS2_blast0.tsv
