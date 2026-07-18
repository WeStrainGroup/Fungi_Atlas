#!/bin/bash
# Purpose: Batch-process single-end atlas sub-projects with the ITS-aware MycoGAP workflow.


PARENT_DIR=""
OUTPUT_DIR=""
THREADS=16

for DIR in $(find $PARENT_DIR -mindepth 1 -maxdepth 1 -type d); do
    PRJ=$(basename $DIR)

    MARKER=$(echo $PRJ | sed 's/.*_//')

    echo "Launching MycoGAP for $PRJ with marker $MARKER ..."

    mkdir -p $OUTPUT_DIR/$PRJ

    nohup mycogap \
        --project e_$PRJ \
        --input ${DIR}/ \
        --seq_type SE \
        --pattern_f .fastq.gz \
        --pattern_r .fastq.gz \
        --marker $MARKER \
        --filter_depth 10000 \
        --thread $THREADS \
        --output $OUTPUT_DIR/$PRJ \
        > $OUTPUT_DIR/$PRJ/e_${PRJ}_log.txt 2>&1 &
done
