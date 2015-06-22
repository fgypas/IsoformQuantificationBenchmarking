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
read_depths <- c(30)
outfolder_plots <- file.path(root, "plots")
prefix <- "sim"
suffix <- "bias_correction.accuracy.trx"
col_select <- c(2,4)
col_names <- c("All transcripts, no bias correction", "All transcripts, bias correction", "Expressed transcripts, no bias correction", "Expressed transcripts, bias correction")
row_select_bias <- c("CEM_bias", "eXpress_bias", "IsoEM_bias", "RSEM_bias", "Sailfish_bias")
row_select_no_bias <- c("CEM_no_bias", "eXpress_no_bias", "IsoEM_no_bias", "RSEM_no_bias", "Sailfish_no_bias")
xlab <- expression("Accuracy (r"[s]*")")

#################
### Functions ###
#################
#<--- Rename methods --->#
rename_methods <- function(names) {
        names <- gsub("_no_bias", "", names)
        names <- gsub("_bias", "", names)
        names <- gsub("$", " (   )", names)
        return(names)
}

############
### Main ###
############
#<--- ITERATE OVER READ DEPTHS --->#
for (read_depth in read_depths) {

        #<--- BUILD MATRIX --->#
        df <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".trx.spearman.tab")), sep="\t")
        df <- df[, col_select]

        #<--- BUILD OUTPUT FILE PREFIX --->#
        outfile <- paste(prefix, read_depth, suffix, sep=".")

        #<--- FILTER METHODS/CONDITIONS CONDITIONS & RENAME --->#
        df_bias <- df[row_select_bias,]
        df_no_bias <- df[row_select_no_bias,]
        rownames(df_bias) <- rename_methods(rownames(df_bias))
        rownames(df_no_bias) <- rename_methods(rownames(df_no_bias))
        
        #<--- MERGE --->#
        df <- merge(df_no_bias, df_bias, by=0, all=TRUE, sort=FALSE)
        rownames(df) <- df$Row.names
        df <- df[,-1]
	df <- df[,c(1,3,2,4)]
        colnames(df) <- col_names

        #<--- PLOT --->#
        dump <- plot_accuracy_per_category(data=df, prefix=file.path(outfolder_plots, outfile), xlab=xlab, colors=c(rep("#DC322f",2), rep("#268BD2", 2)), chars=c(rep(1:2, 2)))

}
