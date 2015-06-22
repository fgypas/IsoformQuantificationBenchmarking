#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 21-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

# rSeq_generate_mapping.py
# rSeq_fix_SAM.pl
# rSeq_calculate_isoform_abundance.py
# samtools
# python 
# perl 

# Shell options
set -e
set -u
set -o pipefail

########################
###  INPUT FILES     ###
########################

transcript_sequences=$1
bam=$2
tmp=$3/rSeq
estimates_dir=$4

mkdir --parents $tmp

########################
###  REFORMAT INPUT  ###
########################

alignments_bam=$tmp/$(basename $bam .bam).ns
# Sort alignment file
samtools sort -n $bam $alignments_bam
# BAM to SAM
sam=$tmp/$(basename $alignments_bam).sam
samtools view -h $alignments_bam.bam > $sam

python rSeq_generate_mapping.py $transcript_sequences $tmp/rSEQ.fa
grep '>' $tmp/rSEQ.fa | awk '{key=substr($0, 19, 15); value=substr($0, 2); dict=key"\t"value; print dict}' > $tmp/rSEQ_dictionary.tab
perl rSeq_fix_SAM.pl $tmp/rSEQ_dictionary.tab $sam > $tmp/rSEQ_fixed.sam

########################
###  RUN TOOL        ###
########################

rseq \
expression_analysis \
$tmp/rSEQ.fa \
$tmp/rSEQ_fixed.sam 

########################
###  ESTIMATES       ###
########################

python rSeq_calculate_isoform_abundance.py --in $tmp/rSEQ_fixed.sam.exp.xls --out $estimates_dir/rSeq.estimates

echo 'rSeq calculation finished !!!'
