#!/bin/bash
# Purpose: Infer continent-level SparCC correlations with FastSpar and estimate p-values from 1,000 bootstraps.


# Parameter input
cd /data/wangxinyu/ITS_Public/Project/Analysis/fastspar
input="/data/wangxinyu/ITS_Public/Project/Analysis/data/location_specific/Asia_otu_gen_p0.01_filterdp.csv"    # Input CSV file path
prj="Asia_gen_0.01"   # Output prefix, e.g., myproject
threads=16
jobs=$((78/threads))

echo $input
echo $prj
echo $threads
echo $jobs

mkdir -p "$prj"
cd "$prj"

start_time=$(date +%s)
echo "[$(date)] FastSpar pipeline started."

# Convert CSV to TSV and add "#OTU ID"
echo "▶ Converting CSV to TSV..."
awk -F',' 'BEGIN {OFS="\t"}
NR==1 { $1 = "#OTU ID"; for (i=1; i<=NF; i++) gsub(/^"|"$/, "", $i); print; next }
      { for (i=1; i<=NF; i++) gsub(/^"|"$/, "", $i); print }
' "$input" > "${prj}_otu.tsv"

# Compute correlation and covariance matrices
echo "▶ Calculating correlation and covariance matrix..."
fastspar --otu_table "${prj}_otu.tsv" \
         --correlation "${prj}_cor.tsv" \
         --covariance "${prj}_cov.tsv" \
         --threads "$threads"

# Bootstrap resampling of OTU table
echo "▶ Bootstrapping OTU table..."
mkdir -p "${prj}_bootstrap_otu"
fastspar_bootstrap --otu_table "${prj}_otu.tsv" \
                   --number 1000 \
                   --prefix "${prj}_bootstrap_otu/boot"

# Compute correlation (cor) and covariance (cov) matrices for these bootstrap tables using GNU parallel
echo "▶ Running fastspar on bootstrap replicates..."
mkdir -p "${prj}_bootstrap_cor"
parallel --jobs $jobs fastspar \
  --otu_table {} \
  --correlation "${prj}_bootstrap_cor/cor_{/}" \
  --covariance "${prj}_bootstrap_cor/cov_{/}" \
  --threads "$threads" \
  ::: "${prj}_bootstrap_otu/"*

# Estimate p-values based on bootstrap correlations
echo "▶ Calculating p-values..."
fastspar_pvalues \
  --otu_table "${prj}_otu.tsv" \
  --correlation "${prj}_cor.tsv" \
  --prefix "${prj}_bootstrap_cor/cor_" \
  --permutations 1000 \
  --outfile "${prj}_pval.tsv"

echo "Done. Output files:"
echo "- Correlation: ${prj}_cor.tsv"
echo "- Covariance: ${prj}_cov.tsv"
echo "- P-values:   ${prj}_pval.tsv"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "[$(date)] FastSpar pipeline completed."
echo "Total runtime: $((runtime / 60)) min $((runtime % 60)) sec"
