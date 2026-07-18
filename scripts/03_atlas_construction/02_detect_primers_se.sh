#!/bin/bash
# Purpose: Batch-detect primer and amplicon structure in single-end projects with APHunter.


PARENT_DIR=""
OUTPUT_DIR=""
THREADS=16

for DIR in $(find $PARENT_DIR -mindepth 1 -maxdepth 1 -type d); do
    PRJ=$(basename $DIR)
    echo "Launching APHunter for $PRJ"
    mkdir -p "$OUTPUT_DIR/$PRJ"
    nohup aphunter -i "${DIR}/" -o "$OUTPUT_DIR/$PRJ" -s .fastq.gz -t $THREADS 2>&1 &
done
