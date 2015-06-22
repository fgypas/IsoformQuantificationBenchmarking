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

infolder <- file.path(root, "analyzed_data")
read_depths <- c(1, 3, 10, 30, 100)
col_select <- 4
select <- c("BitSeq", "CEM", "Cufflinks", "eXpress", "IsoEM", "MMSEQ", "RSEM", "rSeq", "Sailfish", "Scripture", "TIGAR2", "Counting_Transcript", "Counting_Union_exon")
prefix_values <- file.path(root, "plots", "sim.gene.metrics.scaled_values.accuracy_vs_read_depth.expressed")
prefix_ranks <- file.path(root, "plots", "sim.gene.metrics.ranks.accuracy_vs_read_depth.expressed")
palette <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200)[30:170]
col_order <- c(1, 4, 7, 10, 13, 2, 5, 8, 11, 14, 3, 6, 9, 12, 15)
number_color <- "#363636"
number_cex <- 0.9
modified_corrplot <- file.path(root, "scripts", "plotting_functions.corrplot_modified.R")

###############
## Functions ##
###############
## Rename methods
rename_rows <- function(names) {
        names <- gsub("Counting_Union_exon", "Count: 'Union exon'", names)
        names <- gsub("Counting_Transcript", "Count: 'Transcript'", names)
        names <- gsub("_no_priors", " ('priors' to 0)", names)
        return(names)
}
## Rename columns
rename_columns <- function(names) {
	names <- gsub("sim_", "", names)
	names <- gsub(".", " million reads, ", names, fixed=TRUE)
	return(names)
}
## Plot heatmap
plot_heatmap <- function(values.heat, values.num, prefix, width=7, height=7, fontsize.png=16, col.heat=palette, col.num=number_color, cex.num=number_cex) {
	pdf(paste(prefix, "pdf", sep="."), width=width, height=height)
	corrplot(values.heat, is.corr=FALSE, method="shade", col=col.heat, tl.col="black", cl.pos="n")
	corrplot(values.num, is.corr=FALSE, method="number", col=col.num, tl.pos="n", cl.pos="n", number.cex=cex.num, bg="transparent", add=TRUE)
	dev.off()	
	png(paste(prefix, "png", sep="."), width=width, height=height, units="in", res=90, pointsize=fontsize.png)
	corrplot(values.heat, is.corr=FALSE, method="shade", col=col.heat, tl.col="black", cl.pos="n")
	corrplot(values.num, is.corr=FALSE, method="number", col=col.num, tl.pos="n", cl.pos="n", number.cex=cex.num, bg="transparent", add=TRUE)
	dev.off()
	svg(paste(prefix, "svg", sep="."), width=width, height=height)
	corrplot(values.heat, is.corr=FALSE, method="shade", col=col.heat, tl.col="black", cl.pos="n")
	corrplot(values.num, is.corr=FALSE, method="number", col=col.num, tl.pos="n", cl.pos="n", number.cex=cex.num, bg="transparent", add=TRUE)
	dev.off()	
} 

##########
## Main ##
##########
## COMPILE DATA
values_ls <- lapply(read_depths, function(read_depth) {
        #<--- LOAD METRICS --->#
        spearman <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".gene.spearman.tab")), sep="\t")[ , col_select, drop=FALSE]
        pearson <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".gene.pearson.tab")), sep="\t")[ , col_select, drop=FALSE]
        rmse <- read.table(file.path(infolder, paste0("reference_comparison_Ground_truth.sim_", read_depth, ".gene.rmse.tab")), sep="\t")[ , col_select, drop=FALSE]
	#<--- MERGE --->#
	df <- merge(spearman, pearson, by=0, all=TRUE)
	df <- merge(df, rmse, by.x=1, by.y=0, all=TRUE)
	rownames(df) <- df[ , 1]
	df <- df[select, -1]
	colnames(df) <- c("Spearman", "Pearson", "RMSE")
	#<--- RETURN --->#
	return(df)
})

## SET NAMES
names(values_ls) <- paste("sim", read_depths, sep="_")

## SCALE VALUES BETWEEN 0 (WORST) & 1 (BEST)
values_adj_ls <- lapply(values_ls, function(experiment) {
	experiment$Spearman <- experiment$Spearman - min(experiment$Spearman)
	experiment$Spearman <- experiment$Spearman / max(experiment$Spearman)
	experiment$Pearson <- experiment$Pearson - min(experiment$Pearson)
	experiment$Pearson <- experiment$Pearson / max(experiment$Pearson)
	experiment$RMSE <- (experiment$RMSE - max(experiment$RMSE)) * -1
	experiment$RMSE <- experiment$RMSE / max(experiment$RMSE)
	return(experiment)
})

## RANK VALUES
ranks_ls <- lapply(values_adj_ls, function(experiment) {
	as.data.frame(apply(experiment, 2, rank))
})

## TO MATRIX
values <- do.call(cbind, values_ls)
values <- as.matrix(values[, col_order])
values_adj <- do.call(cbind, values_adj_ls)
values_adj <- as.matrix(values_adj[, col_order])
ranks <- do.call(cbind, ranks_ls)
ranks <- as.matrix(ranks[, col_order])

## RENAME ROWS & COLUMNS
colnames(values) <- colnames(values_adj) <- colnames(ranks) <- rename_columns(colnames(values))
rownames(values) <- rownames(values_adj) <- rownames(ranks) <- rename_rows(rownames(values))

## LOAD PLOTTING PACKAGE
library(corrplot)

## LOAD MODIFIED CORRPLOT FUNCTION
source(modified_corrplot)

## PLOT HEATMAP: ADJUSTED VALUES
plot_heatmap(values_adj, values, prefix_values)

## PLOT HEATMAP: RANKS
plot_heatmap(ranks, values, prefix_ranks)
