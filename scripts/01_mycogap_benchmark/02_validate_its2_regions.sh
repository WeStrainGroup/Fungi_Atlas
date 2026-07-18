#!/bin/bash
# Purpose: Validate candidate amplicons and extract fungal ITS2 regions with ITSx.


ITSx \
-i refseq_fun_obipcr_its2.fa \
-o itsx/refseq_fun_obipcr_its2 \
-t F \
--detailed_results T \
--table T \
--cpu 64
