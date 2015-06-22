#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 21-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

module load Java/1.7.0_21
# samtools

# Shell options
set -e
set -u
set -o pipefail

########################
###  INPUT FILES     ###
########################

gtf=$1
bam=$2
tmp=$3/IsoEM
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

alignments_bam=$tmp/$(basename $bam .bam).cs
# Sort alignment file
samtools sort -n $bam $alignments_bam
# BAM to SAM
sam=$tmp/$(basename $alignments_bam .bam).sam
samtools view -h $alignments_bam.bam > $sam

########################
###  RUN TOOL        ###
########################

isoem \
-G $gtf \
-m $mean \
-d $sd \
$sam \
-o $tmp/IsoEM \
--directional

########################
###  ESTIMATES       ###
########################

cp $tmp/IsoEM.iso_estimates $estimates_dir/IsoEM.estimates
cp $tmp/IsoEM.gene_estimates $estimates_dir/IsoEM.estimates.gene

echo 'IsoEM finished !!!'
