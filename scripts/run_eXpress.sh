#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 21-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

module load eXpress/1.5.1-goolf-1.4.10
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
tmp=$3/eXpress
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

express \
$transcript_sequences \
$alignments_bam.bam \
--output-dir $tmp \
--no-update-check \
--f-stranded \
--frag-len-mean $mean \
--frag-len-stddev $sd

########################
###  ESTIMATES       ###
########################

cut -f 2,11 $tmp/results.xprs | tail -n +2 > $estimates_dir/eXpress.estimates

echo 'eXpress finished !!!'
