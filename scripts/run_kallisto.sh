#!/bin/bash

#########################################################
### Foivos Gypas, Biozentrum, University of Basel     ###
### foivos.gypas@unibas.ch                            ###
### 22-MAY-2015                                       ###
#########################################################

########################
###  PRE-REQUISITES  ###
########################

module load Sailfish/0.6.3-goolf-1.4.10-beta

# Shell options
set -e
set -u
set -o pipefail

########################
###  INPUT FILES     ###
########################

transcript_sequences=$1
reads=$2
tmp=$3/Sailfish
estimates_dir=$4
threads=$5

mkdir --parents $tmp

########################
###  RUN TOOL        ###
########################

# Sailfish index
sailfish index \
-t $transcript_sequences \
-o $tmp/Sailfish.idx \
-k 20 \
-p $threads

# Sailfish quantification
sailfish quant \
-i $tmp/Sailfish.idx \
-l T=SE:S=S \
-r $reads \
-o $tmp/Sailfish.quantification \
-p $threads

########################
###  ESTIMATES       ###
########################

cut -f 1,4 $tmp/Sailfish.quantification/quant.sf | tail -n +6 > $estimates_dir/Sailfish.estimates

echo 'Sailfish finished !!!'
