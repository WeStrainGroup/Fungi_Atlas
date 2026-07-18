#!/bin/bash
# Purpose: Generate ITS-aware MycoGAP reference profiles for the WeP1, HMP, and CHGM cohorts.



nohup mycogap \
--project WeP1 \
--input /data/wangxinyu/Fungi_Atlas/Data/WeGut_raw/WeP1/dada/raw \
--pattern_f .1.fq.gz \
--pattern_r .2.fq.gz \
--marker ITS2 \
--filter_depth 10000 \
--thread 8 \
--output /data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/WeP1 \
> /data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/WeP1/log.txt 2>&1 &


nohup mycogap \
--project HMP \
--input /data/wangxinyu/Fungi_Atlas/Data/Public_raw/PE/PRJNA356769b_ITS2 \
--pattern_f .1.fastq.gz \
--pattern_r .2.fastq.gz \
--marker ITS2 \
--filter_depth 10000 \
--thread 8 \
--output /data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/HMP \
> /data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/HMP/log.txt 2>&1 &


nohup mycogap \
--project CHGM \
--input /data/wangxinyu/Fungi_Atlas/Data/Public_raw/PE/PRJCA010668_ITS1 \
--pattern_f .1.fastq.gz \
--pattern_r .2.fastq.gz \
--marker ITS1 \
--filter_depth 10000 \
--thread 8 \
--output /data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/CHGM \
>  /data/wangxinyu/Fungi_Atlas/Analysis/Q_ITS_Error/Annotation/mycogap/CHGM/log.txt 2>&1 &
