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
read_depths <- c(30)
outfolder_data <- file.path(root, "tmp", "plots")
outfolder_plots <- file.path(root, "plots")
prefix <- "sim"
suffix <- "accuracy_per_bin.exons_per_transcript.trx"
col_select <- "^Exons_per_transcript"
row_select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2")
col_names <- c("1 (788)", "2 (1,585)", "3 (1,359)", "4 (1,203)", "5-6 (1,747)", "7-10 (1,385)", "11-107 (1,435)")
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

############
### Main ###
############
#<--- ITERATE OVER READ DEPTHS --->#
for (read_depth in read_depths) {

	#<--- BUILD MATRIX --->#
        df <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".trx.spearman.tab")), sep="\t")
	df <- df[,grep(col_select, colnames(df))]

	#<--- BUILD OUTPUT FILE PREFIX --->#
	outfile <- paste(prefix, read_depth, suffix, sep=".")

        #<--- WRITE OUTPUT TABLE --->#
        write.table(df, file.path(outfolder_data, paste(outfile, "tab", sep=".")), quote=FALSE, sep="\t")

	#<--- SET COLUMN NAMES FOR LEGEND --->#
        colnames(df) <- col_names

	#<--- FILTER METHODS/CONDITIONS CONDITIONS & RENAME --->#
	df <- df[row_select,]
	rownames(df) <- rename_methods(rownames(df))

	#<--- PLOT --->#
	dump <- plot_accuracy_per_subset(data=df, prefix=file.path(outfolder_plots, outfile), xlab=xlab)
}
