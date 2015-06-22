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

#################
### Functions ###
#################
#<--- Rename methods ---># --->#
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
infile_trx <- "replicate_agreement.hsa_no_filter.trx.spearman.tab"
infile_gene <- "replicate_agreement.hsa_no_filter.gene.spearman.tab"
outprefix <- "hsa.trx_gene.replicate_agreement"
select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "Counting_Transcript", "Counting_Union_exon")

#<--- LOAD TRX DATA --->#
trx <- read.table(file.path(infolder, infile_trx))[ , 1, drop=FALSE]
colnames(trx) <- "Transcripts"

#<--- LOAD GENE DATA --->#
gene <- read.table(file.path(infolder, infile_gene))[ , 1, drop=FALSE]
colnames(gene) <- "Genes"

#<--- MERGE --->#
corr_df <- merge(trx, gene, by=0, all=TRUE)
rownames(corr_df) <- corr_df[,1]
corr_df <- corr_df[,-1]

#<--- SELECT METHODS/CONDITIONS --->#
corr_df <- corr_df[select,]

#<--- RENAME METHODS/CONDITIONS --->#
rownames(corr_df) <- rename_methods(rownames(corr_df))

#<--- PLOT --->#
dump <- plot_agreement_per_category(corr_df, file.path(outfolder_plots, outprefix), min_x=0, max_x=1)

#############
### Mouse ###
#############
#<--- SPECIFIC PARAMETERS --->#
infile_trx <- "replicate_agreement.mmu_no_filter.trx.spearman.tab"
infile_gene <- "replicate_agreement.mmu_no_filter.gene.spearman.tab"
outprefix <- "mmu.trx_gene.replicate_agreement"
select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "Counting_Transcript", "Counting_Union_exon")

#<--- LOAD TRX DATA --->#
trx <- read.table(file.path(infolder, infile_trx))[ , 1, drop=FALSE]
colnames(trx) <- "Transcripts"

#<--- LOAD GENE DATA --->#
gene <- read.table(file.path(infolder, infile_gene))[ , 1, drop=FALSE]
colnames(gene) <- "Genes"

#<--- MERGE --->#
corr_df <- merge(trx, gene, by=0, all=TRUE)
rownames(corr_df) <- corr_df[,1]
corr_df <- corr_df[,-1]

#<--- SELECT METHODS/CONDITIONS --->#
corr_df <- corr_df[select,]

#<--- RENAME METHODS/CONDITIONS --->#
rownames(corr_df) <- rename_methods(rownames(corr_df))

#<--- PLOT --->#
dump <- plot_agreement_per_category(corr_df, file.path(outfolder_plots, outprefix), min_x=0, max_x=1)
