#!/bin/bash

SOURCE_DIR=$1
TARGET_DIR=$2

cd $SOURCE_DIR

find . -name 'accepted_hits.bam'
inputfiles=$(find . -name 'accepted_hits.bam' | xargs echo)
readlink -m $inputfiles

featureCounts -a hg19.gtf -o rnaseq.txt $inputfiles
