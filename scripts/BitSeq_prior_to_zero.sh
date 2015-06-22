#!/bin/bash

## Alexander Kanitz
## 21-JAN-2015

in_folder=$1
out_file=$2
tmp_folder=$3

# Set temp folder to current folder if not specified
if [[ "$tmp_folder" = "" ]]; then
	tmp_folder=$PWD
fi

# Check for and save to variable ID file path (*.trx)
if [ -f $in_folder/*.trx_ids ]; then
        id_file=$in_folder/*.trx_ids
else
        >&2 echo "There is either no or more than one '*.mean.rpkm' file in the input folder. Execution aborted!"; exit
fi

# Check for and save to variable RPKM file path (*.mean.rpkm)
if [ -f $in_folder/*.mean.rpkm ]; then
	rpkm_file=$in_folder/*.mean.rpkm
else
	>&2 echo "There is either no or more than one '*.mean.rpkm' file in the input folder. Execution aborted!"; exit
fi

# Check for and save to variable estimated count file path (*.m_alphas)
if [ -f $in_folder/*.m_alphas ]; then
        count_file=$in_folder/*.m_alphas
else
        >&2 echo "There is either no or more than one '*.mean.rpkm' file in the input folder. Execution aborted!"; exit
fi

# Make random-name folder in temp folder
tmp_folder=`mktemp -d $tmp_folder/tmp.XXXXXX`

# Copy transcript IDs to temp folder
cp $id_file $tmp_folder/ids

# Extract RPKM column to file in temp folder
tail -n +4 $rpkm_file | cut --delim " " -f1 > $tmp_folder/rpkm

# Extract estimated counts column to file in temp folder
tail -n +8 $count_file | cut --delim " " -f2 > $tmp_folder/counts

# Paste IDs, RPKMs and counts
paste $tmp_folder/ids $tmp_folder/rpkm $tmp_folder/counts > $tmp_folder/combined

# Set priors to zeros and write to output file
awk -v OFS="\t" '{if ($3 == 1) print $1, 0; else print $1, $2}' $tmp_folder/combined > $out_file

# Remove temporary files
rm -rf $tmp_folder
