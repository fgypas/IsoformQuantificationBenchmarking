#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 20-MAY-2015                                       ###
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

gtf=$1
bam=$2
tmp=$3/Cufflinks
estimates_dir=$4
threads=$5

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
samtools sort $bam $alignments_bam

########################
###  RUN TOOL        ###
########################

module load Cufflinks/2.1.1-goolf-1.4.10

time cufflinks \
--output-dir $tmp \
--GTF $gtf \
-p $threads \
--library-type fr-secondstrand \
$alignments_bam.bam \
--multi-read-correct \
--frag-len-mean $mean \
--frag-len-std-dev $sd


########################
###  ESTIMATES       ###
########################

tail -n +2 $tmp/isoforms.fpkm_tracking | grep -v "^ENSG" | grep -P "\tOK$" | cut -f1,10 > $estimates_dir/Cufflinks.estimates
tail -n +2 $tmp/genes.fpkm_tracking | cut -f1,10 > $estimates_dir/Cufflinks.estimates.gene

echo 'Cufflinks finished !!!' 
