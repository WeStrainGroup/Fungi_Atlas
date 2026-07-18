#!/bin/bash
# Purpose: Sample 50 reference sequences into each of 50 reproducible mock communities.


#INPUT_FILE="/home/wangxy/ITS_benchmark/refseq_ITS/refseq_fun_obipcr_its2_derep.fa"
INPUT_FILE="/home/wangxy/ITS_benchmark/data/refseq_fun_obipcr_its2_v97_filter.fa"
OUTPUT_DIR="/home/wangxy/ITS_benchmark/data/sampling"


for i in {1..50}
do
    printf -v j "%02d" $i

    echo "doing $j ..."
    seqkit sample2 -n 50 -s $i ${INPUT_FILE} > "$OUTPUT_DIR/mock_${j}.fa"
done

echo "all done!"
