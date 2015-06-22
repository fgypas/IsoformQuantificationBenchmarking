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
infolder <- file.path(root, "tmp", "plots")
read_depths <- c(30)
outfolder_data <- file.path(root, "tmp", "plots")
outfolder_plots <- file.path(root, "plots")
prefix <- "sim"
suffix <- "accuracy_per_bin.standard_deviations"
row_select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2")
col_names <- c("Transcript length", "GC content", "Exons per transcript", "Transcripts per gene")
xlab <- expression(sigma ~ "accuracy (r"[s]*")")

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
	sd <- list()
	sd$`Transcript_length` <- read.table(file.path(infolder, paste("sim", read_depth, "accuracy_per_bin", "transcript_length", "trx", "tab", sep=".")), sep="\t")
	sd$`GC_content` <- read.table(file.path(infolder, paste("sim", read_depth, "accuracy_per_bin", "gc_content", "trx", "tab", sep=".")), sep="\t")
	sd$`Exons_per_transcript` <- read.table(file.path(infolder, paste("sim", read_depth, "accuracy_per_bin", "exons_per_transcript", "trx", "tab", sep=".")), sep="\t")
	sd$`Transcripts_per_gene` <- read.table(file.path(infolder, paste("sim", read_depth, "accuracy_per_bin", "transcripts_per_gene", "gene", "tab", sep=".")), sep="\t")
	sd <- lapply(sd, apply, 1, sd)
	sd$`Transcripts_per_gene` <- sd$`Transcripts_per_gene`[! names(sd$`Transcripts_per_gene`) %in% c("Counting_Union_exon", "Counting_Transcript")]
	df <- do.call(cbind, sd)

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
	dump <- plot_accuracy_per_subset(data=df, prefix=file.path(outfolder_plots, outfile), xlab=xlab, max_x=0.3)

}
