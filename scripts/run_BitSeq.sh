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

# Transcript sequence files
trx_fa=$1
# Alignments
bam=$2
# Temporary directory
tmp=$3/BitSeq
# Outfile prefix
out_prefix=$tmp/BitSeq
# Estimates
estimates_dir=$4
# Threads
threads=$5

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

# BitSeq parseAlignment
parseAlignment $alignments_bam.bam --outFile ${out_prefix}.prob --trSeqFile $trx_fa --trInfoFile ${out_prefix}.trx --uniform -P $threads --veryVerbose &>${out_prefix}.parseAlignment.log
# BitSeq estimateVBExpression
estimateVBExpression ${out_prefix}.prob --outPrefix ${out_prefix} --outType RPKM --trInfoFile ${out_prefix}.trx --samples 1000 --seed 1 -P $threads --veryVerbose &>${out_prefix}.estimateVBExpression.log
# BitSeq getVariance
getVariance ${out_prefix}.VBrpkm --outFile ${out_prefix}.mean.rpkm --verbose &>${out_prefix}.getVariance.log
# Get transcript IDs
cut --delim " " -f2 ${out_prefix}.trx | tail -n +2 > ${out_prefix}.trx_ids
# Merge transcript IDs and RPKM values
cut --delim " " -f1 ${out_prefix}.mean.rpkm | tail -n +4 | paste ${out_prefix}.trx_ids - > ${estimates_dir}/BitSeq.estimates

echo 'BitSeq finished !!!' 
