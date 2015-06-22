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
#<--- Rename methods --->#
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
infile_trx <- "replicate_agreement.hsa.trx.spearman.tab"
infile_pas <- "replicate_agreement.hsa.pas.spearman.tab"
infile_gene <- "replicate_agreement.hsa.gene.spearman.tab"
outprefix_all <- "hsa.trx_pas_gene.replicate_agreement.all_features"
outprefix_expressed <- "hsa.trx_pas_gene.replicate_agreement.expressed_features"
select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "A-seq", "Counting_Transcript", "Counting_Union_exon")

#<--- LOAD TRX DATA --->#
all_trx <- read.table(file.path(infolder, infile_trx))[ , 4, drop=FALSE]
expressed_trx <- read.table(file.path(infolder, infile_trx))[ , 5, drop=FALSE]

#<--- LOAD PAS DATA --->#
all_pas <- read.table(file.path(infolder, infile_pas))[ , 2, drop=FALSE]
expressed_pas <- read.table(file.path(infolder, infile_pas))[ , 4, drop=FALSE]

#<--- LOAD GENE DATA --->#
all_gene <- read.table(file.path(infolder, infile_gene))[ , 4, drop=FALSE]
expressed_gene <- read.table(file.path(infolder, infile_gene))[ , 6, drop=FALSE]

#<--- MERGE FEATURES CORRESPONDING TO ALL POLY(A) SITES --->#
mt_all <- merge(all_trx, all_pas, by=0, all=TRUE)
mt_all <- merge(mt_all, all_gene, by.x=1, by.y=0, all=TRUE)
dimnames(mt_all) <- list(mt_all$Row.names, c("Row.names", "Transcripts", "3'-end processing sites", "Genes"))

#<--- MERGE FEATURES CORRESPONDING TO EXPRESSED POLY(A) SITES --->#
mt_expressed <- merge(expressed_trx, expressed_pas, by=0, all=TRUE)
mt_expressed <- merge(mt_expressed, expressed_gene, by.x=1, by.y=0, all=TRUE)
dimnames(mt_expressed) <- list(mt_expressed$Row.names, c("Row.names", "Transcripts", "3'-end processing sites", "Genes"))

#<--- SELECT & RENAME METHODS/CONDITIONS --->#
mt_all <- mt_all[select,-1]
rownames(mt_all) <- rename_methods(rownames(mt_all))
mt_expressed <- mt_expressed[select,-1]
rownames(mt_expressed) <- rename_methods(rownames(mt_expressed))

#<--- PLOT --->#
dump <- plot_agreement_per_category(mt_all, file.path(outfolder_plots, outprefix_all), min_x=0, max_x=1)
dump <- plot_agreement_per_category(mt_expressed, file.path(outfolder_plots, outprefix_expressed), min_x=0, max_x=1)

#############
### Mouse ###
#############
#<--- SPECIFIC PARAMETERS --->#
infile_trx <- "replicate_agreement.mmu.trx.spearman.tab"
infile_pas <- "replicate_agreement.mmu.pas.spearman.tab"
infile_gene <- "replicate_agreement.mmu.gene.spearman.tab"
outprefix_all <- "mmu.trx_pas_gene.replicate_agreement.all_features"
outprefix_expressed <- "mmu.trx_pas_gene.replicate_agreement.expressed_features"
select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "A-seq", "Counting_Transcript", "Counting_Union_exon")

#<--- LOAD TRX DATA --->#
all_trx <- read.table(file.path(infolder, infile_trx))[ , 4, drop=FALSE]
expressed_trx <- read.table(file.path(infolder, infile_trx))[ , 5, drop=FALSE]

#<--- LOAD PAS DATA --->#
all_pas <- read.table(file.path(infolder, infile_pas))[ , 2, drop=FALSE]
expressed_pas <- read.table(file.path(infolder, infile_pas))[ , 4, drop=FALSE]

#<--- LOAD GENE DATA --->#
all_gene <- read.table(file.path(infolder, infile_gene))[ , 4, drop=FALSE]
expressed_gene <- read.table(file.path(infolder, infile_gene))[ , 6, drop=FALSE]

#<--- MERGE FEATURES CORRESPONDING TO ALL POLY(A) SITES --->#
mt_all <- merge(all_trx, all_pas, by=0, all=TRUE)
mt_all <- merge(mt_all, all_gene, by.x=1, by.y=0, all=TRUE)
dimnames(mt_all) <- list(mt_all$Row.names, c("Row.names", "Transcripts", "3'-end processing sites", "Genes"))

#<--- MERGE FEATURES CORRESPONDING TO EXPRESSED POLY(A) SITES --->#
mt_expressed <- merge(expressed_trx, expressed_pas, by=0, all=TRUE)
mt_expressed <- merge(mt_expressed, expressed_gene, by.x=1, by.y=0, all=TRUE)
dimnames(mt_expressed) <- list(mt_expressed$Row.names, c("Row.names", "Transcripts", "3'-end processing sites", "Genes"))

#<--- SELECT & RENAME METHODS/CONDITIONS --->#
mt_all <- mt_all[select,-1]
rownames(mt_all) <- rename_methods(rownames(mt_all))
mt_expressed <- mt_expressed[select,-1]
rownames(mt_expressed) <- rename_methods(rownames(mt_expressed))

#<--- PLOT --->#
dump <- plot_agreement_per_category(mt_all, file.path(outfolder_plots, outprefix_all), min_x=0, max_x=1)
dump <- plot_agreement_per_category(mt_expressed, file.path(outfolder_plots, outprefix_expressed), min_x=0, max_x=1)
