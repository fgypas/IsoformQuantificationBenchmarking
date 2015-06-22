#!/bin/bash

########################################################################
### Foivos Gypas						     ###
### foivos.gypas@unibas.ch				             ###
### 26-MAY-2015   						     ###
### Biozentrum, University of Basel                                  ###
########################################################################

# AK: This is currently only a wrapper for estimating the expression from a single sample. The real wrapper 4C would call this wrapper 9 times...
# AK: Guess what's missing here still is command-line parsing (together with DOC/usage & some basic validation) to get root, transcriptome, gtf, alignments etc. when called from the superwrapper.
# AK: I have modified the calls a bit. I hope they still work... I explicitly specified the 'bash' shell, since I'm not sure we can guarantee that it would work with whatever shell people have defined as 'sh' (e.g. if you specified 'set -e' in any of the individual method wrappers, things would fail already for the 'zsh' etc). The requirement for 'bash' shell is explicitly stated in '0.README', so that should be fine like that. I have also wrapped the variables in quotes since it is not impossible that someone might have a space in the path to the root directory.

# FG: Should I write a different script for Sailfish bias? I think not sice the sailfish calculates the results with and without bias at the same time and reports 2 files. 

set -e 

#######################################
## PARAMETERS SIMULATIONS 1 MILLION  ##
#######################################

root=$1
transcriptome=$2
gtf=$3
alignment_transcriptome=$4
alignment_genome=$5
tmp=$5
estimates=$6
threads=$7
trx_gene_lookup_table=$8
reads=$9 # AK: Guess we still have to copy the FASTA/FASTQ files to $root/scripts/seq_libs

# RUN BitSeq
bash "$root/scripts/run_BitSeq.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates" "$threads"

# RUN CEM 
bash "$root/scripts/run_CEM.sh" "$gtf" "$alignment_genome" "$tmp" "$estimates"

# RUN Cufflinks
bash "$root/scripts/run_Cufflinks.sh" "$gtf" "$alignment_genome" "$tmp" "$estimates" "$threads"

# RUN eXpress
bash "$root/scripts/run_eXpress.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates"

# RUN IsoEM
bash "$root/scripts/run_IsoEM.sh" "$gtf" "$alignment_genome" "$tmp" "$estimates"
 
# RUN MMSEQ
bash "$root/scripts/run_MMSEQ.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates"

# RUN RSEM
bash "$root/scripts/run_RSEM.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates" "$threads" "$trx_gene_lookup_table"

# RUN rSeq
bash "$root/scripts/run_rSEQ.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates"

# RUN Sailfish
bash "$root/scripts/run_Sailfish.sh" "$transcriptome" "$reads" "$tmp" "$estimates" "$threads"

# RUN Scripture
bash "$root/scripts/run_Scripture.sh" "$gtf" "$alignment_genome" "$tmp" "$estimates"

# RUN TIGAR2 
bash "$root/scripts/run_TIGAR2.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates"

# RUN CEM bias
bash "$root/scripts/run_CEM_bias.sh" "$gtf" "$alignment_genome" "$tmp" "$estimatesOC"

# RUN IsoEM bias
bash "$root/scripts/run_IsoEM_bias.sh" "$gtf" "$alignment_genome" "$tmp" "$estimates"

# RUN RSEM bias
bash "$root/scripts/run_RSEM_bias.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates" "$threads" "$trx_gene_lookup_table"

# RUN eXpress no bias
bash "$root/scripts/run_eXpress_no_bias.sh" "$transcriptome" "$alignment_transcriptome" "$tmp" "$estimates"
