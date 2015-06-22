#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 22-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

module load Java/1.7.0_21

# Shell options
set -e
set -u
set -o pipefail

########################
###  INPUT FILES     ###
########################

transcript_sequences=$1
bam=$2
tmp=$3/TIGAR2
estimates_dir=$4

mkdir --parents $tmp

########################
###  REFORMAT INPUT  ###
########################

alignments_bam=$tmp/$(basename $bam .bam).ns
# Sort alignment file
samtools sort -n $bam $alignments_bam

########################
###  RUN TOOL        ###
########################

time java -Xms16g -Xmx16g -jar ~/bin/Tigar2.jar \
$transcript_sequences \
$alignments_bam.bam \
--alpha_zero 0.1 \
$tmp/TIGAR2

########################
###  ESTIMATES       ###
########################

cut -f 1,4 $tmp/TIGAR2 | tail -n +2 > $estimates_dir/TIGAR2.estimates

echo 'TIGAR2 finished !!!'
