#!/bin/bash

## Alexander Kanitz
## 16-FEB-2015

## In the output of the MMSEQ method for the inference of gene/isoform abundance from RNA-seq data (file *.mmseq), resets the estimates for all features whose standard deviation of the posteriors is bigger than the maximum posterior standard deviation of all features that have at least one unique hit to "NA".

## Used procedure ("unique hits")
tail -n +3 $1 | awk -F"\t" -v OFS="\t" '{ if ($8 == 0) $2 = "NA" ; else exp($2); print $1, $2 }'

## Alternative procedures ("observed")
#tail -n +3 $1 | awk -F"\t" -v OFS="\t" '{ if ($13 == 0) $2 = "NA" ; print $1, $2 }'
