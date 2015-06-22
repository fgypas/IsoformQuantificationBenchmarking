#!/usr/bin/Rscript

############
## Header ##
############
# Created: Jan 02, 2015
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
subsets <- file.path(root, "resources", "feature_subsets.sim.R")
sim <- file.path(root, "estimates", "sim_1.transcripts.Ground_truth.reference")
lookup <- file.path(root, "resources", "gencode.v19.annotation.ENS_compatible.trx_gene_lookup_table")
prefix <- file.path(root, "plots", "sim.expression_cdf")
colors <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#000000")

###############
## Load data ##
###############
load(subsets)
expr <- read.table(sim, col.names=c("trx_id", "value"), stringsAsFactors=FALSE, colClasses=c("character", "numeric"))
trx_2_gene <- read.table(lookup, col.names=c("trx_id", "gene_id"), stringsAsFactors=FALSE, colClasses=c("character", "character"))

###################
## Load packages ##
###################
if ( suppressWarnings(suppressPackageStartupMessages(require("KernSmooth"))) == FALSE ) { stop("Package 'KernSmooth' required!\nExecution aborted.") }

#####################
## Data processing ##
#####################
# Get/add gene expression
expr <- merge(trx_2_gene, expr)
expr_gene <- aggregate(value~gene_id, expr, sum)
expr <- setNames(c(expr$value, expr_gene$value), c(expr$trx_id, expr_gene$gene_id))
# Set pseudocount and take log
expr[expr == 0] <- 1/32
expr <- log2(expr)
# Get IDs per bin
metrics <- list()
metrics$`Transcript length` <- subsets$trx[9:12]
metrics$`GC content` <- subsets$trx[20:23]
metrics$`Exons per transcript` <- subsets$trx[13:19]
metrics$`Transcripts per gene` <- subsets$gene[9:13]
# Rename bins
metrics <- lapply(metrics, function(metric) {
	names(metric) <- sapply(strsplit(names(metric), "--col--_"), "[", 2)
	names(metric) <- gsub("_to_", "-", names(metric))
	names(metric) <- gsub("_", " ", names(metric))
	names(metric) <- gsub("--op--", "(", names(metric))
	names(metric) <- gsub("--cl--", " transcripts)", names(metric))
	names(metric) <- gsub("--com--", ",", names(metric))
	names(metric) <- gsub(" nts ", " nt ", names(metric))
	return(metric)
})
names(metrics$`GC content`) <- gsub(" (", "% (", names(metrics$`GC content`), fixed=TRUE)
names(metrics$`Transcripts per gene`) <- gsub("transcripts)", "genes)", names(metrics$`Transcripts per gene`), fixed=TRUE) 
# Get expression values * cumulative distribution function per bin
expr_bin <- lapply(metrics, lapply, function(bin) expr[names(expr) %in% bin])
ecdf_bin <- lapply(expr_bin, lapply, ecdf)

##########
## Plot ##
##########
# Loop over different metrics
for (metric in names(ecdf_bin)) {
	# Get minimum and maximum expression values for current metric
	max <- round(max(unlist(expr_bin[[metric]])), digits=1)
	min <- round(min(unlist(expr_bin[[metric]])), digits=1)	
	#<--- PLOT PDF --->#
	# Open plot device and re-size margins
	pdf(paste(prefix, gsub(" ", "_", tolower(metric)), "pdf", sep="."))
	par(mar=c(5.1,5.1,4.1,1.1))
	# Plot frame and abline
	plot(0, type="n", xlim=c(min,max), ylim=c(0,1), main=metric, xlab="Log2 TPM", ylab="Cumulative fraction", cex.axis=1.25, cex.lab=1.25, cex.main=1.5)
	abline(0.5, 0, lty=2)
	# Plot legend
	legend_x <- max - max(strwidth(names(ecdf_bin[[metric]]), units="user")) - 0.12 * abs(max-min)
	legend_y <- max(strheight(names(ecdf_bin[[metric]]), units="user")) * (length(ecdf_bin[[metric]]) * 1.65 + 1) + 0.02
	legend(legend_x, legend_y, legend=names(ecdf_bin[[metric]]), col=colors[1:length(ecdf_bin[[metric]])], lwd=1.4, box.lwd=0, box.col="white", bg="white")
	# Plot data
	for (bin in 1:length(ecdf_bin[[metric]])) {
		plot(ecdf_bin[[metric]][[bin]], add=TRUE, lwd=1.25, col=colors[[bin]], do.points=FALSE)
	}
	# Close plot device
	graphics.off()
	#<--- PLOT SVG --->#
	# Open plot device and re-size margins
	svg(paste(prefix, gsub(" ", "_", tolower(metric)), "svg", sep="."))
	par(mar=c(5.1,5.1,4.1,1.1))
	# Plot frame and abline
	plot(0, type="n", xlim=c(min,max), ylim=c(0,1), main=metric, xlab="Log2 TPM", ylab="Cumulative fraction", cex.axis=1.25, cex.lab=1.25, cex.main=1.5)
	abline(0.5, 0, lty=2)
	# Plot legend
	legend_x <- max - max(strwidth(names(ecdf_bin[[metric]]), units="user")) - 0.12 * abs(max-min)
	legend_y <- max(strheight(names(ecdf_bin[[metric]]), units="user")) * (length(ecdf_bin[[metric]]) * 1.65 + 1) + 0.02
	legend(legend_x, legend_y, legend=names(ecdf_bin[[metric]]), col=colors[1:length(ecdf_bin[[metric]])], lwd=1.4, box.lwd=0, box.col="white", bg="white")
	# Plot data
	for (bin in 1:length(ecdf_bin[[metric]])) {
		plot(ecdf_bin[[metric]][[bin]], add=TRUE, lwd=1.25, col=colors[[bin]], do.points=FALSE)
	}
	# Close plot device
	graphics.off()
	#<--- PLOT PNG --->#
	# Open plot device and re-size margins
	png(paste(prefix, gsub(" ", "_", tolower(metric)), "png", sep="."), width=480, height=480, pointsize=15)
	par(mar=c(5.1,5.1,4.1,1.1))
	# Plot frame and abline
	plot(0, type="n", xlim=c(min,max), ylim=c(0,1), main=metric, xlab="Log2 TPM", ylab="Cumulative fraction", cex.axis=1.25, cex.lab=1.25, cex.main=1.5)
	abline(0.5, 0, lty=2)
	# Plot legend
	legend_x <- max - max(strwidth(names(ecdf_bin[[metric]]), units="user", font=15)) - 0.15 * abs(max-min)
	legend_y <- max(strheight(names(ecdf_bin[[metric]]), units="user", font=15)) * (length(ecdf_bin[[metric]]) * 1.65 + 1) + 0.02
	legend(legend_x, legend_y, legend=names(ecdf_bin[[metric]]), col=colors[1:length(ecdf_bin[[metric]])], lwd=1.4, box.lwd=0, box.col="white", bg="white")
	# Plot data
	for (bin in 1:length(ecdf_bin[[metric]])) {
		plot(ecdf_bin[[metric]][[bin]], add=TRUE, lwd=1.25, col=colors[[bin]], do.points=FALSE)
	}
	# Close plot device
	graphics.off()
}
