#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 20-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

# gtf_transcripts_to_bed12.R
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
tmp=$3/CEM
estimates_dir=$4

mkdir --parents $tmp

########################
###  REFORMAT INPUT  ###
########################

# Generate bed from gtf
cem_bed=$tmp/$(basename $gtf)'.bed'
Rscript gtf_transcripts_to_bed12.R -i $gtf -o $cem_bed
# Remove first line
tail -n +2 $cem_bed > $cem_bed.trimmed
# Sort by chromosome and position
sort -k 1,1 -k 2,2n $cem_bed.trimmed > $cem_bed.sorted

# Sort alignment file
samtools view -H $bam > $tmp/$(basename $bam .bam).head
samtools view $bam | sort -T $tmp -k 3,3 -k 4,4n | cat "$tmp/$(basename $bam .bam).head" - | samtools view -bS - > $tmp/CEM.bam

########################
###  RUN TOOL        ###
########################

# Run CEM
runcem.py \
--annotation $cem_bed.sorted \
--forceref \
--prefix $tmp/CEM \
$tmp/CEM.bam

########################
###  ESTIMATES       ###
########################

cut -f 2,9 $tmp/CEM*.pred | tail -n +3 > $tmp/CEM.tmp
head -n -1 $tmp/CEM.tmp > $estimates_dir/CEM.estimates
rm $tmp/CEM.tmp
