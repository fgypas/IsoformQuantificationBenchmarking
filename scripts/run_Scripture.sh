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
module load SAMtools/0.1.19-goolf-1.4.10
module load BEDTools/2.18.1-goolf-1.4.10

# Shell options
set -e
set -u
set -o pipefail

########################
###  INPUT FILES     ###
########################

gtf=$1
bam=$2
tmp=$3/Scripture
estimates_dir=$4

mkdir --parents $tmp

########################
###  REFORMAT INPUT  ###
########################

alignments_bam=$tmp/$(basename $bam .bam).cs
# Sort alignment file
samtools sort $bam $alignments_bam

scripture_bed=$tmp/Scripture.bed
scripture_alignment=$tmp/Scripture.list
echo $alignments_bam.bam > $scripture_alignment
samtools index $alignments_bam.bam 
Rscript gtf_transcripts_to_bed12.R -i $gtf -o $scripture_bed

########################
###  RUN TOOL        ###
########################

java -Xmx16g -jar ~/bin/ScriptureScorer.jar \
-alignments $scripture_alignment \
-annotations $scripture_bed \
-out $tmp/Scripture.exp \
-minMappingQuality 5 -strand first -singleEnd

########################
###  ESTIMATES       ###
########################

tail -n +2 $tmp/Scripture.exp.rpkm.txt > $estimates_dir/Scripture.estimates

echo 'Scripture finished !!!'
