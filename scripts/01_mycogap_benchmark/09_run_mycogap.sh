#!/bin/bash
# Purpose: Process all simulated mock communities with MycoGAP in paired-end ITS2 mode.


    nohup time mycogap \
        --project mock \
        --input /home/wangxy/ITS_benchmark/data/fq \
        --pattern_f _R1.fastq.gz \
        --pattern_r _R2.fastq.gz \
        --marker ITS2 \
        --filter_depth 10000 \
        --thread 16 \
        --output /home/wangxy/ITS_benchmark/process/mycogap \
        > /home/wangxy/ITS_benchmark/process/mycogap/log_mycogap.txt 2>&1 &
