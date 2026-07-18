#!/bin/bash
# Purpose: Align ASVs within each selected SH and infer bootstrap-supported phylogenetic trees.


BASE="/data/wangxinyu/ITS_Public/Project/Analysis/genetic_distance"
CLUSTER_DIR="$BASE/vsearch/clusters"
TREE_DIR="$BASE/tree"
CSV="$BASE/clustered/SHs_selected.csv"

# Read SH column from CSV, skipping header
tail -n +2 "$CSV" | cut -d',' -f2 | tr -d '"' | while read -r SH; do
  # Generate file paths
  FASTA="$CLUSTER_DIR/${SH}.fasta"
  ALIGN_DIR="$TREE_DIR/${SH}"
  ALIGN_FILE="$ALIGN_DIR/${SH}_align.fasta"

  # Create output directory
  mkdir -p "$ALIGN_DIR"

  # Run MUSCLE alignment
  echo "Aligning $FASTA ..."
  muscle --align "$FASTA" -output "$ALIGN_FILE" -threads 64

  # Build phylogenetic tree
  echo "Running IQ-TREE on $ALIGN_FILE ..."
  nohup iqtree -s "$ALIGN_FILE" -m MFP -B 1000 -bnni -T AUTO > "$ALIGN_DIR/iqtree.log" 2>&1 &
done
