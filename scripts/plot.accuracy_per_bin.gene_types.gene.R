#!/usr/bin/Rscript

############
## Header ##
############
# Created: Jan 04, 2015
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
plot_function_source <- file.path(root, "scripts", "plotting_functions.accuracy_per_subset.R")
source(plot_function_source)
infolder <- file.path(root, "analyzed_data")
outfolder_data <- file.path(root, "tmp", "plots")
outfolder_plots <- file.path(root, "plots")
col_select <- "^Gene_biotype"
xlab <- expression("Accuracy (r"[s]*")")

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
infile1 <- "reference_comparison_A-seq-2.hsa_1.gene.spearman.tab"
infile2 <- "reference_comparison_A-seq-2.hsa_2.gene.spearman.tab"
prefix <- "hsa.accuracy_per_bin.gene_types.gene"
row_select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "Counting_Transcript", "Counting_Union_exon")
col_names <- c("Protein-coding (12,513)", "lincRNA (735)", "Antisense (739)")

#<--- BUILD MATRIX --->#
df1 <- read.table(file.path(infolder, infile1), sep="\t")
df2 <- read.table(file.path(infolder, infile2), sep="\t")
df1 <- df1[rownames(df1) %in% rownames(df2), ]
df2 <- df2[rownames(df2) %in% rownames(df1), ]
df1 <- df1[match(rownames(df1), rownames(df2)), ]
df <- (df1 + df2) / 2
df <- df[,grep(col_select, colnames(df))]

#<--- WRITE OUTPUT TABLE --->#
write.table(df, file.path(outfolder_data, paste(prefix, "tab", sep=".")), quote=FALSE, sep="\t")

#<--- SET COLUMN NAMES FOR LEGEND --->#
colnames(df) <- col_names

#<--- FILTER METHODS/CONDITIONS & RENAME --->#
df <- df[row_select,]
rownames(df) <- rename_methods(rownames(df))

#<--- PLOT --->#
dump <- plot_accuracy_per_subset(data=df, prefix=file.path(outfolder_plots, prefix), xlab=xlab, colors=c("#DC322f", "#268BD2", "#009914"))

#############
### Mouse ###
#############
#<--- SPECIFIC PARAMETERS --->#
infile1 <- "reference_comparison_A-seq-2.mmu_1.gene.spearman.tab"
infile2 <- "reference_comparison_A-seq-2.mmu_2.gene.spearman.tab"
prefix <- "mmu.accuracy_per_bin.gene_types.gene"
row_select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "Counting_Transcript", "Counting_Union_exon")
col_names <- c("Protein-coding (9,639)", "lincRNA (91)", "Antisense (55)")

#<--- BUILD MATRIX --->#
df1 <- read.table(file.path(infolder, infile1), sep="\t")
df2 <- read.table(file.path(infolder, infile2), sep="\t")
df1 <- df1[rownames(df1) %in% rownames(df2), ]
df2 <- df2[rownames(df2) %in% rownames(df1), ]
df1 <- df1[match(rownames(df1), rownames(df2)), ]
df <- (df1 + df2) / 2
df <- df[,grep(col_select, colnames(df))]

#<--- WRITE OUTPUT TABLE --->#
write.table(df, file.path(outfolder_data, paste(prefix, "tab", sep=".")), quote=FALSE, sep="\t")

#<--- SET COLUMN NAMES FOR LEGEND --->#
colnames(df) <- col_names

#<--- FILTER METHODS/CONDITIONS & RENAME --->#
df <- df[row_select,]
rownames(df) <- rename_methods(rownames(df))

#<--- PLOT --->#
dump <- plot_accuracy_per_subset(data=df, prefix=file.path(outfolder_plots, prefix), xlab=xlab, colors=c("#DC322f", "#268BD2", "#009914"))
