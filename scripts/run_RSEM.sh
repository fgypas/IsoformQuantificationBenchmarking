#!/bin/bash

#########################################################
### Foivos Gypas, Alexander Kanitz		      ###
### Biozentrum, University of Basel     	      ###
### foivos.gypas@unibas.ch                            ###
### 22-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

# samtools
# python

# Shell options
set -e
set -u
set -o pipefail

########################
###  INPUT FILES     ###
########################

transcript_sequences=$1
bam=$2
tmp=$3/RSEM
estimates_dir=$4
threads=$5
look_up_table=$6

mkdir --parents $tmp

#######################
## PREPARE ANNOTATION #
######################

cut -f1 $look_up_table > $tmp/annotation.trx_ids
cut -f2 $look_up_table > $tmp/annotation.gene_ids 
paste $tmp/annotation.gene_ids $tmp/annotation.trx_ids > $tmp/annotation.gene_trx_lookup_table
rsem_index_files_prefix=$tmp/$(basename $transcript_sequences .fa)

rsem-prepare-reference --no-polyA --transcript-to-gene-map $tmp/annotation.gene_trx_lookup_table $transcript_sequences $rsem_index_files_prefix > $tmp/index.log

########################
###  REFORMAT INPUT  ###
########################

alignments_bam=$tmp/$(basename $bam .bam).ns
# Sort alignment file
samtools sort -n $bam $alignments_bam

## Variables for filtering reads by CIGAR string
awk="/bin/awk"
sam_trx_rsem="$tmp/trx_RSEM.sam"

# Filter read file: remove reads with INDELs and other prohibited characters in the CIGAR string
samtools view -h $alignments_bam.bam | $awk '{if ($6 !~ /[DHINPS]/ ) print $0}' > $sam_trx_rsem

## Variables for recalculating read length distribution (mean & SD)
sam_trx_rsem_mean="$tmp/mean"
sam_trx_rsem_sd="$tmp/sd"

# Calculate read length distribution
python sam_read_length_stats_no_pipe.py --sam $sam_trx_rsem --mean $sam_trx_rsem_mean --sd $sam_trx_rsem_sd --multimappers

## Variables for "rsem-calculate-expression" wrapper script
rsem_seed_length=15
rsem_fragment_length_mean=`cat $sam_trx_rsem_mean`
rsem_fragment_length_sd=`cat $sam_trx_rsem_sd`
rsem_output_files_prefix="$tmp/RSEM_results"

########################
###  RUN TOOL        ###
########################

# Build command for "rsem-calculate-expression"
rsem_command="rsem_calculate_expression --sam --strand-specific --no-qualities --seed-length $rsem_seed_length --fragment-length-mean $rsem_fragment_length_mean --fragment-length-sd $rsem_fragment_length_sd -p $threads --time $sam_trx_rsem $rsem_index_files_prefix $rsem_output_files_prefix"

$rsem_command

########################
###  ESTIMATES       ###
########################

cut -f 1,7 $tmp/RSEM_results*.isoforms.results | tail -n +2 > $estimates_dir/RSEM.estimates
cut -f 1,7 $tmp/RSEM_results*.genes.results | tail -n +2 > $estimates_dir/RSEM.estimates.gene

echo 'RSEM finished !!!'
