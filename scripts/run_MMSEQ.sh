#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 21-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

# samtools

# Shell options
set -e
set -u
set -o pipefail

########################
###  INPUT FILES     ###
########################

transcript_sequences=$1
bam=$2
tmp=$3/MMSEQ
estimates_dir=$4

mkdir --parents $tmp

########################
### CALCULATE STATS  ###
########################

samtools view $bam | python sam_read_length_stats.py --mean $tmp/mean --sd $tmp/sd --multimappers
mean=$(cat $tmp/mean)
sd=$(cat $tmp/sd)

########################
###  REFORMAT INPUT  ###
########################

alignments_bam=$tmp/$(basename $bam .bam).ns
# Sort alignment file
samtools sort -n $bam $alignments_bam

########################
###  RUN TOOL        ###
########################

# Bam2Hits
bam2hits-linux \
$transcript_sequences \
$alignments_bam.bam > $tmp/FINAL.hits

# MMSEQ
mmseq-linux \
$tmp/FINAL.hits \
$tmp/FINAL.hits.bam2hits

########################
###  ESTIMATES       ###
########################
cut -f1,2 $tmp/*.bam2hits.mmseq | tail -n +3 > $tmp/MMSEQ.logmu
Rscript tab_exponentiate_column_tab.R --input $tmp/MMSEQ.logmu --output $estimates_dir/MMSEQ.estimates

cut -f1,2 $tmp/*.bam2hits.gene.mmseq | tail -n +3 > $tmp/MMSEQ_gene.logmu
Rscript tab_exponentiate_column_tab.R --input $tmp/MMSEQ_gene.logmu --output $estimates_dir/MMSEQ.estimates.gene

echo 'MMSEQ finished !!!'
