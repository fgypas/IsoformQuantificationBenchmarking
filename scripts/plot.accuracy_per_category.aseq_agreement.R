#!/usr/bin/Rscript

############
## Header ##
############
# Created: Jan 05, 2015
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
colnames <- c("Transcripts, replicate 1", "Transcripts, replicate 2", "3'-end processing sites, replicate 1", "3'-end processing sites, replicate 2", "Genes, replicate 1", "Genes, replicate 2")
select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "Counting_Transcript", "Counting_Union_exon")

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
        names <- gsub("Counting_Union_exon", "Count: 'Union exon'", names)
        names <- gsub("Counting_Transcript", "Count: 'Transcript'", names)
        names <- gsub("_no_priors", " ('priors' to 0)", names)
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
outprefix <- "hsa.trx_pas_gene.aseq_agreement"

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
colnames(mt) <- colnames

#<--- SELECT METHODS/CONDITIONS --->#
mt <- mt[select,]

#<--- RENAME METHODS/CONDITIONS --->#
rownames(mt) <- rename_methods(rownames(mt))

#<--- PLOT --->#
dump <- plot_accuracy_per_category(mt, file.path(outfolder_plots, outprefix), min_x=0, max_x=1, colors=c(rep("#DC322f", 2), rep("#268BD2", 2), rep("#009914", 2)), chars=rep(1:2, 3))

#############
### Mouse ###
#############
#<--- SPECIFIC PARAMETERS --->#
infile_pas_1 <- "reference_comparison_A-seq-2.mmu_1.pas.spearman.tab"
infile_pas_2 <- "reference_comparison_A-seq-2.mmu_2.pas.spearman.tab"
infile_gene_1 <- "reference_comparison_A-seq-2.mmu_1.gene.spearman.tab"
infile_gene_2 <- "reference_comparison_A-seq-2.mmu_2.gene.spearman.tab"
outprefix <- "mmu.trx_pas_gene.aseq_agreement"

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
colnames(mt) <- colnames

#<--- SELECT METHODS/CONDITIONS --->#
mt <- mt[select,]

#<--- RENAME METHODS/CONDITIONS --->#
rownames(mt) <- rename_methods(rownames(mt))

#<--- PLOT --->#
dump <- plot_accuracy_per_category(mt, file.path(outfolder_plots, outprefix), min_x=0, max_x=1, colors=c(rep("#DC322f", 2), rep("#268BD2", 2), rep("#009914", 2)), chars=rep(1:2, 3))
