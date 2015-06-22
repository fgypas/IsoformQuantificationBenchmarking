#!/usr/bin/Rscript

############
## Header ##
############
# Created: Feb 21, 2015
# Author: Alexander Kanitz
# Company: Zavolan Group, Biozentrum, University of Basel

###################
## CLI arguments ##
###################
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]

################
## Parameters ##
################
plot_function_source <- file.path(root, "scripts", "plotting_functions.false_pos_neg_vs_read_depth.R")
source(plot_function_source)
infolder <- file.path(root, "analyzed_data")
read_depths <- c(1, 3, 10, 30, 100)
col_select <- 5
select <- c("BitSeq", "BitSeq_no_priors", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "MMSEQ_no_priors", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2")
prefix_false_pos <- file.path(root, "plots", "sim.trx.false_pos_vs_read_depth.low_expressed")
prefix_false_neg <- file.path(root, "plots", "sim.trx.false_neg_vs_read_depth.low_expressed")
prefix_true_pos <- file.path(root, "plots", "sim.trx.true_pos_vs_read_depth.low_expressed")
prefix_true_neg <- file.path(root, "plots", "sim.trx.true_neg_vs_read_depth.low_expressed")

###############
## Functions ##
###############
## Rename methods
rename_methods <- function(names) {
        names <- gsub("Counting_Union_exon", "Count: 'Union exon'", names)
        names <- gsub("Counting_Transcript", "Count: 'Transcript'", names)
        names <- gsub("_no_priors", " ('priors' to 0)", names)
        return(names)
}
## Merge multiple data frames from list
merge_df_from_ls <- function(ls, by=0, all=TRUE) {
        mdf <- ls[[1]]
        for (df in 2:length(ls)) {
                mdf <- merge(mdf, ls[[df]], by=by, all=all)
                rownames(mdf) <- mdf$Row.names
                mdf <- mdf[,-1]
        }
        return(mdf)
}

##########
## Main ##
##########
## COMPILE DATA
ls_false_pos <- lapply(read_depths, function(read_depth) {
        df <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".trx.false_pos_rate.tab")), sep="\t")[ , col_select, drop=FALSE]
	colnames(df) <- paste("sim", read_depth, sep="_")
	return(df)
})
ls_false_neg <- lapply(read_depths, function(read_depth) {
        df <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".trx.false_neg_rate.tab")), sep="\t")[ , col_select, drop=FALSE]
        colnames(df) <- paste("sim", read_depth, sep="_")
        return(df)
})
ls_true_pos <- lapply(read_depths, function(read_depth) {
        df <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".trx.true_pos_rate.tab")), sep="\t")[ , col_select, drop=FALSE]
        colnames(df) <- paste("sim", read_depth, sep="_")
        return(df)
})
ls_true_neg <- lapply(read_depths, function(read_depth) {
        df <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".trx.true_neg_rate.tab")), sep="\t")[ , col_select, drop=FALSE]
        colnames(df) <- paste("sim", read_depth, sep="_")
        return(df)
})

## MERGE DATA FRAMES
df_false_pos <- merge_df_from_ls(ls_false_pos)
df_false_neg <- merge_df_from_ls(ls_false_neg)
df_true_pos <- merge_df_from_ls(ls_true_pos)
df_true_neg <- merge_df_from_ls(ls_true_neg)

## SELECT METHODS
df_false_pos <- df_false_pos[select, ]
df_false_neg <- df_false_neg[select, ]
df_true_pos <- df_true_pos[select, ]
df_true_neg <- df_true_neg[select, ]

## RENAME METHODS
rownames(df_false_pos) <- rename_methods(rownames(df_false_pos))
rownames(df_false_neg) <- rename_methods(rownames(df_false_neg))
rownames(df_true_pos) <- rename_methods(rownames(df_true_pos))
rownames(df_true_neg) <- rename_methods(rownames(df_true_neg))

## PLOT
plot_false_neg_depth(df_false_neg, prefix_false_neg, 0.0, 1.0)
plot_true_pos_depth(df_true_pos, prefix_true_pos, 0.0, 1.0)
