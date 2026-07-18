#!/bin/bash
# Purpose: Build/download reference databases and BLAST atlas ITS and sampled 16S sequences against UNITE, NCBI, and SILVA.

conda install bioconda::blast
conda install conda-forge::ncbi-datasets-cli
conda install bioconda::taxonkit
conda install -c bioconda -c conda-forge -c defaults csvtk # Specify channel order


# Build UNITE BLAST database
makeblastdb -in /data/wangxinyu/ITS_Public/Project/Analysis/ref/unite_blast_db/sh_general_release_dynamic_19.02.2025.fasta -dbtype nucl -out /data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/unite_blast_db/unite_blast_db

# Download NCBI database
cd /data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/ncbi_ITS_blast_db
update_blastdb.pl --showall
update_blastdb.pl --decompress ITS_RefSeq_Fungi

# Download NCBI taxonomy dump for taxonkit annotation
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

mkdir -p $HOME/.taxonkit
cp names.dmp nodes.dmp delnodes.dmp merged.dmp $HOME/.taxonkit


# Dereplicate fungal ASVs
cd /data/wangxinyu/ITS_Public/Project/Analysis/data/ASV
vsearch --derep_fulllength refseq_filterdp.fasta --output refseq_filterdp_derep.fasta --sizeout --uc refseq_filterdp_derep.uc


# Run blast
# NCBI
export BLASTDB=/data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/ncbi_ITS_blast_db
blastn -query /data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/refseq_filterdp_derep.fasta
       -db ITS_RefSeq_Fungi
       -out /data/wangxinyu/ITS_Public/Project/Analysis/blast/fungi/refseq_filterdp_derep_ncbi.tsv
       -outfmt '6 qseqid sseqid staxids sscinames sblastnames pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore'
       -max_target_seqs 1
       -num_threads 64


# Further annotate taxon names
cd /data/wangxinyu/ITS_Public/Project/Analysis/blast/fungi
paste refseq_filterdp_derep_ncbi.tsv <(cut -f3 refseq_filterdp_derep_ncbi.tsv | sed 's/;.*//' | taxonkit reformat -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}") | csvtk add-header -t -n "qseqid,sseqid,staxids,sscinames,sblastnames,pident,length,qcovs,mismatch,gapopen,qstart,qend,sstart,send,sstrand,evalue,bitscore,taxids,kingdom,phylum,class,order,family,genus,species" > refseq_filterdp_derep_ncbi_final.tsv


# UNITE
export BLASTDB=/data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/unite_blast_db
blastn -query /data/wangxinyu/ITS_Public/Project/Analysis/data/ASV/refseq_filterdp_derep.fasta
       -db unite_blast_db
       -out /data/wangxinyu/ITS_Public/Project/Analysis/blast/fungi/refseq_filterdp_derep_unite.tsv
       -outfmt '6 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore'
       -max_target_seqs 1
       -num_threads 64

# Add header
cd /data/wangxinyu/ITS_Public/Project/Analysis/blast/fungi
(echo -e "qseqid\tsseqid\tpident\tlength\tqcovs\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue\tbitscore"; cat refseq_filterdp_derep_unite.tsv) > refseq_filterdp_derep_unite_final.tsv



# Use 16S dataset from Cell 2025 as reference
# Download NCBI 16S database
update_blastdb.pl --decompress 16S_ribosomal_RNA

# Build SILVA database
makeblastdb -in /data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/silva_blast_db/SILVA_138.2_SSURef_NR99_tax_silva.fasta -dbtype nucl -out /data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/silva_blast_db/silva_blast_db

# Sample 5% of filtered ASVs after filtering
cd /data/wangxinyu/ITS_Public/Project/Analysis/blast/cell2025/analysis
seqkit sample -p 0.05 cell2025_asv_filter2.fasta > cell2025_asv_filter2_0.05.fasta

# Run blast
# For large datasets, use the run_blast script for parallel computation to save time
# NCBI
export BLASTDB=/data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/ncbi_16S_blast_db
nohup blastn \
-query cell2025_asv_filter2_0.05.fasta \
-db 16S_ribosomal_RNA \
-out cell2025_asv_filter2_0.05_ncbi.tsv \
-outfmt '6 qseqid sseqid staxids sscinames sblastnames pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore' \
-max_target_seqs 1 \
-num_threads 64 > blast_ncbi_log.txt 2>&1 &


# Taxonomic annotation + header
paste cell2025_asv_filter2_0.05_ncbi.tsv <(cut -f3 cell2025_asv_filter2_0.05_ncbi.tsv | sed 's/;.*//' | taxonkit reformat -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}") | csvtk add-header -t -n "qseqid,sseqid,staxids,sscinames,sblastnames,pident,length,qcovs,mismatch,gapopen,qstart,qend,sstart,send,sstrand,evalue,bitscore,taxids,kingdom,phylum,class,order,family,genus,species" > cell2025_asv_filter2_0.05_ncbi_final.tsv


# SILVA
export BLASTDB=/data/wangxinyu/ITS_Public/Project/Analysis/blast/ref/silva_blast_db
nohup blastn \
-query cell2025_asv_filter2_0.05.fasta \
-db silva_blast_db \
-out cell2025_asv_filter2_0.05_silva.tsv \
-outfmt '6 qseqid sseqid pident length qcovs mismatch gapopen qstart qend sstart send sstrand evalue bitscore' \
-max_target_seqs 1 \
-num_threads 64 > blast_silva_log.txt 2>&1 &


# Add header
(echo -e "qseqid\tsseqid\tpident\tlength\tqcovs\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tsstrand\tevalue\tbitscore"; cat cell2025_asv_filter2_0.05_silva.tsv) > cell2025_asv_filter2_0.05_silva_final.tsv
