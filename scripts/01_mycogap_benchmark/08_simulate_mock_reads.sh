#!/bin/bash
# Purpose: Simulate paired-end MiSeq reads from each mock sequence and abundance profile with InSilicoSeq.




INPUT_DIR="/home/wangxy/ITS_benchmark/data/fq"


mkdir -p "$INPUT_DIR"


for input_file in ${INPUT_DIR}/*.fa; do

    base_name=$(basename "$input_file" _sequence.fa)
    abundance=${INPUT_DIR}/${base_name}_abundance.txt

    echo "doing: $base_name"


    iss generate \
        --genomes "$input_file" \
        --abundance_file "$abundance" \
        --sequence_type amplicon \
        --model miseq \
        --n_reads 200000 \
        --compress \
        --cpus 32 \
        --output "${INPUT_DIR}/${base_name}"

    echo "done: $base_name"
done

rm ${INPUT_DIR}/*.vcf

echo "all done!"
