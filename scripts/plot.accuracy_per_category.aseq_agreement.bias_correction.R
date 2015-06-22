#!/usr/bin/Rscript

############
## Header ##
############
# Created: Feb 20, 2015
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
plot_function_source <- file.path(root, "scripts", "plotting_functions.accuracy_per_category.R")
source(plot_function_source)
infolder <- file.path(root, "analyzed_data")
outfolder_plots <- file.path(root, "plots")
col_select_trx <- 5
col_select_pas <- 2
col_select_gene <- 4
colnames <- c("Transcripts", "3'-end processing sites", "Genes")
row_select_bias <- c("CEM_bias", "eXpress_bias", "IsoEM_bias", "RSEM_bias", "Sailfish_bias")
row_select_no_bias <- c("CEM_no_bias", "eXpress_no_bias", "IsoEM_no_bias", "RSEM_no_bias", "Sailfish_no_bias")
xlab <- expression("Accuracy (r"[s]*")")

#################
### Functions ###
#################
## Merge multiple dataframes from list
merge_df_from_ls <- function(ls, by=0, all=TRUE) {
        mdf <- ls[[1]]
        for (df in 2:length(ls)) {
                mdf <- merge(mdf, ls[[df]], by=by, all=all)
                rownames(mdf) <- mdf$Row.names
                mdf <- mdf[,-1]
        }
        return(mdf)
}
## Rename methods
rename_methods <- function(names) {
        names <- gsub("_no_bias", "", names)
        names <- gsub("_bias", "", names)
	names <- gsub("$", " (   )", names)
        return(names)
}

#############
### Human ###
#############
#<--- SPECIFIC PARAMETERS --->#
infile_pas_1 <- "reference_comparison_A-seq-2.hsa_1.pas.spearman.tab"
infile_pas_2 <- "reference_comparison_A-seq-2.hsa_2.pas.spearman.tab"
infile_gene_1 <- "reference_comparison_A-seq-2.hsa_1.gene.spearman.tab"
infile_gene_2 <- "reference_comparison_A-seq-2.hsa_2.gene.spearman.tab"
outprefix_1 <- "hsa.rep_1.bias_correction.accuracy"
outprefix_2 <- "hsa.rep_2.bias_correction.accuracy"

#<--- INITIALIZE DATA CONTAINER LIST --->#
ls <- list()

#<--- LOAD TRX DATA --->#
ls$trx_1 <- read.table(file.path(infolder, infile_pas_1))[ , col_select_trx, drop=FALSE]
ls$trx_2 <- read.table(file.path(infolder, infile_pas_2))[ , col_select_trx, drop=FALSE]

#<--- LOAD PAS DATA --->#
ls$pas_1 <- read.table(file.path(infolder, infile_pas_1))[ , col_select_pas, drop=FALSE]
ls$pas_2 <- read.table(file.path(infolder, infile_pas_2))[ , col_select_pas, drop=FALSE]

#<--- LOAD GENE DATA --->#
ls$gene_1 <- read.table(file.path(infolder, infile_gene_1))[ , col_select_gene, drop=FALSE]
ls$gene_2 <- read.table(file.path(infolder, infile_gene_2))[ , col_select_gene, drop=FALSE]

#<--- MERGE DATA --->#
mt <- as.matrix(merge_df_from_ls(ls))

#<--- SELECT METHODS/CONDITIONS --->#
mt_bias <- mt[row_select_bias,]
mt_no_bias <- mt[row_select_no_bias,]

#<--- RENAME METHODS/CONDITIONS --->#
rownames(mt_bias) <- rename_methods(rownames(mt_bias))
rownames(mt_no_bias) <- rename_methods(rownames(mt_no_bias))

#<--- MERGE BIAS/NO BIAS ROWS --->#
mt <- merge(mt_no_bias, mt_bias, by=0, all=TRUE, sort=FALSE)
rownames(mt) <- mt[,1]

#<--- SPLIT BY REPLICATE --->#
mt_1 <- mt[ , seq(2, ncol(mt), 2)]
mt_2 <- mt[ , seq(3, ncol(mt), 2)]

#<--- REARRANGE COLUMNS --->#
mt_1 <- mt_1[, c(1,4,2,5,3,6)]
mt_2 <- mt_2[, c(1,4,2,5,3,6)]

#<--- RENAME COLUMNS --->#
colnames(mt_1) <- colnames(mt_2) <- paste0(unlist(lapply(colnames, rep, 2)), rep(c(", no bias correction", ", bias correction"), 3))

#<--- BUILD OUTPUT FILE PREFIXES --->#
outfile_1 <- file.path(outfolder_plots, outprefix_1)
outfile_2 <- file.path(outfolder_plots, outprefix_2)

#<--- PLOT --->#
dump <- plot_accuracy_per_category(data=mt_1, prefix=outfile_1, xlab=xlab, colors=c(rep("#DC322f", 2), rep("#268BD2", 2), rep("#009914", 2)), chars=c(rep(1:2, 3)))
dump <- plot_accuracy_per_category(data=mt_2, prefix=outfile_2, xlab=xlab, colors=c(rep("#DC322f", 2), rep("#268BD2", 2), rep("#009914", 2)), chars=c(rep(1:2, 3)))

#############
### Mouse ###
#############
#<--- SPECIFIC PARAMETERS --->#
infile_pas_1 <- "reference_comparison_A-seq-2.mmu_1.pas.spearman.tab"
infile_pas_2 <- "reference_comparison_A-seq-2.mmu_2.pas.spearman.tab"
infile_gene_1 <- "reference_comparison_A-seq-2.mmu_1.gene.spearman.tab"
infile_gene_2 <- "reference_comparison_A-seq-2.mmu_2.gene.spearman.tab"
outprefix_1 <- "mmu.rep_1.bias_correction.accuracy"
outprefix_2 <- "mmu.rep_2.bias_correction.accuracy"

#<--- INITIALIZE DATA CONTAINER LIST --->#
ls <- list()

#<--- LOAD TRX DATA --->#
ls$trx_1 <- read.table(file.path(infolder, infile_pas_1))[ , col_select_trx, drop=FALSE]
ls$trx_2 <- read.table(file.path(infolder, infile_pas_2))[ , col_select_trx, drop=FALSE]

#<--- LOAD PAS DATA --->#
ls$pas_1 <- read.table(file.path(infolder, infile_pas_1))[ , col_select_pas, drop=FALSE]
ls$pas_2 <- read.table(file.path(infolder, infile_pas_2))[ , col_select_pas, drop=FALSE]

#<--- LOAD GENE DATA --->#
ls$gene_1 <- read.table(file.path(infolder, infile_gene_1))[ , col_select_gene, drop=FALSE]
ls$gene_2 <- read.table(file.path(infolder, infile_gene_2))[ , col_select_gene, drop=FALSE]

#<--- MERGE DATA --->#
mt <- as.matrix(merge_df_from_ls(ls))

#<--- SELECT METHODS/CONDITIONS --->#
mt_bias <- mt[row_select_bias,]
mt_no_bias <- mt[row_select_no_bias,]

#<--- RENAME METHODS/CONDITIONS --->#
rownames(mt_bias) <- rename_methods(rownames(mt_bias))
rownames(mt_no_bias) <- rename_methods(rownames(mt_no_bias))

#<--- MERGE BIAS/NO BIAS ROWS --->#
mt <- merge(mt_no_bias, mt_bias, by=0, all=TRUE, sort=FALSE)
rownames(mt) <- mt[,1]

#<--- SPLIT BY REPLICATE --->#
mt_1 <- mt[ , seq(2, ncol(mt), 2)]
mt_2 <- mt[ , seq(3, ncol(mt), 2)]

#<--- REARRANGE COLUMNS --->#
mt_1 <- mt_1[, c(1,4,2,5,3,6)]
mt_2 <- mt_2[, c(1,4,2,5,3,6)]

#<--- RENAME COLUMNS --->#
colnames(mt_1) <- colnames(mt_2) <- paste0(unlist(lapply(colnames, rep, 2)), rep(c(", no bias correction", ", bias correction"), 3))

#<--- BUILD OUTPUT FILE PREFIXES --->#
outfile_1 <- file.path(outfolder_plots, outprefix_1)
outfile_2 <- file.path(outfolder_plots, outprefix_2)

#<--- PLOT --->#
dump <- plot_accuracy_per_category(data=mt_1, prefix=outfile_1, xlab=xlab, colors=c(rep("#DC322f", 2), rep("#268BD2", 2), rep("#009914", 2)), chars=c(rep(1:2, 3)))
dump <- plot_accuracy_per_category(data=mt_2, prefix=outfile_2, xlab=xlab, colors=c(rep("#DC322f", 2), rep("#268BD2", 2), rep("#009914", 2)), chars=c(rep(1:2, 3)))
