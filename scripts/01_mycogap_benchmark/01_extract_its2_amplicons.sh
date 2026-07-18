#!/bin/bash
# Purpose: Extract candidate ITS2 amplicons from fungal RefSeq sequences with OBITools in-silico PCR.


obipcr \
--forward GCATCGATGAAGAACGCAGC \
--reverse TCCTCCGCTTATTGATATGC \
--max-length 1000 \
--min-length 100 \
--allowed-mismatches 3 \
--delta 0 \
/home/wangxy/ITS_benchmark/data/fungi.ITS.fna \
> fungi.ITS.obipcr.fa
