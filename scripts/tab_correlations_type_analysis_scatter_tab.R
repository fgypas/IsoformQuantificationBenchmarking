#!/usr/bin/Rscript

#=============#
#  HEADER //  #
#=============#
### Created: May 31, 2015
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
#=============#
#  // HEADER  #
#=============#


#==============#
#  OPTIONS //  #
#==============#
#---> LOAD OPTIONS PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Calculate different correlation metrics between expression estimates and generates scatter plots. Two analyses types are supported: (1) omparison of (one or many different) estimates against one or more references, and (2) agreement between replicates.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 31-MAY-2015"
version <- "Version: 1.0.1 (31-MAY-2015)"
requirements <- "Requires: optparse, LSD (if '--plot' is specified), KernSmooth (if '--plot is specified)"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
	make_option("--experiment", action="store", type="character", default=NULL, help="Experiment name (required).", metavar="string"),
	make_option("--feature-type", action="store", type="character", default=NULL, help="Type of the features provided in the reference and sample. One of 'trx' (transcripts), 'pas' (3'-end/polyA processing sites) or 'gene' (genes) (required).", metavar="string"),
	make_option("--feature-subsets-R", action="store", type="character", default=NULL, help="R object, output of 'compile_feature_subsets.[hsa|mmu|sim].R' scripts (required).", metavar="path"),
	make_option("--input-directory", action="store", type="character", default=NULL, help="Directory containing the expression estimates to compare (required). If '--replicate-directory' is not specified, all files matching '--glob-sample' will be compared to all files matching '--glob-reference' with the metrics specified under '--metrics'. One output table will be produced per metric and per reference.", metavar="path"),
	make_option("--reference-directory", action="store", type="character", default=NULL, help="Directory containing the reference estimates. If not specified, this defaults to the argument to '--input-directory'. Only one of '--reference-directory' and '--replicate-directory' may be specified.", metavar="path"),	
	make_option("--replicate-directory", action="store", type="character", default=NULL, help="Directory containing expression estimates computed for a replicate sample. If specified, estimate files matching '--glob-sample' in '--input-directory' will not be compared to references, but rather to the corresponding estimate files in this directory. Filenames of replicate estimates have to be identical. Not compatible with '--reference-directory' and '--glob-reference'.", metavar="path"),
	make_option("--output-directory", action="store", type="character", default=NULL, help="Directory where output metric tables shall be written.", metavar="path"),
	make_option("--output-directory-plots", action="store", type="character", default=NULL, help="Directory where output plots are generated. If not specified, will default to '--input-directory'. Ignored if '--plot' is not specified.", metavar="path"),
	make_option("--glob-sample", action="store", type="character", default="*", help="File selector glob for samples (default='*'). File basenames up to the first '.' are used for as row names in output files.", metavar="glob"),
	make_option("--glob-reference", action="store", type="character", default=NULL, help="File selector glob for references (default='*'). File basenames up to the first '.' are used for generating output filenames. Only one of '--glob-reference' and '--replicate-directory' may be specified.", metavar="glob"),
	make_option("--metrics", action="store", type="character", default="all", help="Metric to use for the comparison of estimates/references. Supported values: 'pearson', 'spearman', 'rmse' (root mean square error), and 'meanratio' (sample over reference). Metrics 'pearson', 'rmse' and 'meanratio' are calculated after setting the pseudocount and conversion of values to log2 space. Multiple metrics may be specified (separated by commas but not spaces). Alternatively, specify 'all' (default) to compute all supported metrics.", metavar="string"),
	make_option("--pseudocount", action="store", type="character", default=1/32, help="Positive number that will be used for replacing any values smaller than this value prior to calculating Pearson correlation coefficients, root mean square errors and plotting. Accepts floating point numbers and fractions (default: '1/32', corresponding to log2 of '-5').", metavar="float"),
	make_option("--true-false-analysis", action="store_true", default=FALSE, help="Analyze true/false positive/negative rates when comparing samples against references. Ignored if '--replicate-directory' is specified."), 
	make_option("--subset-stats", action="store_true", default=FALSE, help="Compute the Kolmogorov-Smirnov test between the sample over reference ratios for all subsets relative to a defined reference subset (see '--subset-stats-references'). Returns one table each for the D statistic of the test and the associated P value."),
	make_option("--subset-stats-references", action="store", type="character", default=NULL, help="Subset to use as references when computing the Kolmogorov-Smirnov test (see '--subset-stats'). Can be either a list index number (integer; default: 1) or an element name. The number(name) will be used to select a reference subset of IDs of the current feature type ('trx', 'pas', 'gene'; see '--feature-type') in the subset container object specified via '--feature-subsets-R'.", metavar="string"),
	make_option("--plot", action="store_true", default=FALSE, help="Generate log2 scatter plots (requires packages 'LSD' and 'KernSmooth')."),
	make_option("--plot-format", action="store", type="character", default="all", help="Output format for scatter plots. Supported formats: 'pdf', 'svg', 'png', 'png_data_only'. Multiple metrics may be specified (separated by commas but not spaces). Alternatively, specify 'all' (default) to generate the plots in all formats. Note that for 'png_data_only' images, only data is plotted (no axes, no labels).", metavar="string"),
	make_option("--plot-min", action="store", type="numeric", default=NULL, help="Minimum limit of plotting range (default: determined by data).", metavar="float"),
	make_option("--plot-max", action="store", type="numeric", default=NULL, help="Maximum limit of plotting range (default: determined by data).", metavar="float"),
	make_option("--help", action="store_true", default=FALSE, help="Show this information and die."),
	make_option("--usage", action="store_true", default=FALSE, dest="help", help="Show this information and die."),
	make_option("--verbose", action="store_true", default=FALSE, help="Print log messages.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --experiment <STRING> --feature-type <STRING> --feature-subsets-R <PATH> --input-directory <PATH> --output-directory <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Re-assign command-line arguments
exp <- opt$`experiment`
type <- opt$`feature-type`
subsetsFile <- opt$`feature-subsets-R`
inDir <- opt$`input-directory`
refDir <- opt$`reference-directory`
repDir <- opt$`replicate-directory`
outDir <- opt$`output-directory`
outDirPlot <- opt$`output-directory-plots`
globSample <- opt$`glob-sample`
globRef <- opt$`glob-reference`
metrics <- opt$`metrics`
trueFalse <- opt$`true-false-analysis`
subsetStats <- opt$`subset-stats`
statsRefs <- opt$`subset-stats-references`
plot <- opt$`plot`
plotFormat <- opt$`plot-format`
plotMin <- opt$`plot-min`
plotMax <- opt$`plot-max`
pseudo <- eval(parse(text=opt$`pseudocount`))
verbose <- opt$`verbose`

## Derive switches
doReps <- ! is.null(repDir)

## Set allowed values
allowedTypes <- c("trx", "pas", "gene")
allowedMetrics <- c("pearson", "spearman", "rmse", "meanratio")
allowedFormats <- c("pdf", "svg", "png", "png_data_only")

## Die if required arguments are missing
if ( is.null(exp) || is.null(type) || is.null(subsetsFile) || is.null(inDir) || is.null(outDir) ) {
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}
## Die if mutually exclusive options are specified
if ( ( ! is.null(refDir) || ! is.null(globRef) ) && doReps ) {
	write("[ERROR] Mutually exclusive options specified! Options '--reference-directory' and '--glob-reference' are not allowed when '--replicate-directory' is specified.\n\n", stderr())
	stop(print_help(opt_parser))
}
## Set defaults if not specified
if ( is.null(refDir) ) refDir <- inDir
if ( is.null(globRef) ) globRef <- "*"
if ( is.null(outDirPlot) ) outDirPlot <- outDir
## Validate allowed values
if ( ! type %in% allowedTypes ) {
	write("[ERROR] Unsupported values specified! Specify one of the following to option '--metrics': 'trx', 'pas', 'gene'.\n\n", stderr())
	stop(print_help(opt_parser))	
}
if ( metrics == "all" ) {
	metrics <- allowedMetrics
} else if ( metrics == "" ) {
	metrics <- metrics
} else {
	metrics <- unlist(strsplit(metrics,","))
	if ( ! all(metrics %in% allowedMetrics) ) {
		write("[ERROR] Unsupported value specified! Specify either 'all' or one or more (separated by comma but no spaces) of the following to option '--metrics': 'pearson', 'spearman', 'rmse'.\n\n", stderr())
		stop(print_help(opt_parser))
	}
}
if ( ! pseudo > 0 ) {
	write("[ERROR] Unsupported values specified! Option '--plot-pseudocount' requires a positive number as its argument.\n\n", stderr())
	stop(print_help(opt_parser))
}
if ( subsetStats ) {
	if ( is.null(statsRefs) ) {
		statsRefs <- 1
	} else {
		int <- suppressWarnings(as.integer(statsRefs))
		statsRefs <- ifelse(is.na(int), statsRefs, int)
	}
}
if ( plotFormat == "all" ) {
	plotFormat <- allowedFormats
} else {
	plotFormat <- unlist(strsplit(plotFormat,","))
	if ( ! all(plotFormat %in% allowedFormats) ) {
		write("[ERROR] Unsupported value specified! Specify either 'all' or one or more (separated by comma but no spaces) of the following to option '--plot-format': 'pdf', 'svg', 'png', 'png_data_only'.\n\n", stderr())
		stop(print_help(opt_parser))
	}
}
#==============#
#  // OPTIONS  #
#==============#


#================#
#  FUNCTIONS //  #
#================#
#---> CALCULATE METRICS FOR REPLICATE AGREEMENTS <---#
calculateReplicateAgreement <- function(samples, replicates, outDir, outDirPlot, experiment, featureType, subsets, subsetStats, statsRefs, metrics, pseudo, plot, format, min, max) {
	# Build output filename prefix
	prefix <- file.path(outDir, paste("replicate_agreement", experiment, featureType, sep="."))
	plotPrefix <- file.path(outDirPlot, paste("replicate_agreement", experiment, featureType, sep="."))
	# Iterate over samples
	results <- mapply(function(sampleName, sample, replicate) {
		# Iterate over subsets
		lapply(names(subsets), function(subset) {
			# Extract sample and replicate values
			m <- merge(sample[sample$id %in% subsets[[subset]], ], replicate[replicate$id %in% subsets[[subset]], ], by=1)
			x <- m[ , 2]	# Data from '--input-directory' on the x axis ("replicate 1")
			y <- m[ , 3]	# Data from '--replicate-directory' on the y axis ("replicate 2")
			# Initialize container for results
			resultsList <- list()
			# Calculate metrics
			if ("spearman" %in% metrics) resultsList$spearman <- cor(x, y, method="spearman")
			x[x < pseudo] <- pseudo
			y[y < pseudo] <- pseudo
			x <- log2(x)
			y <- log2(y)
			if ("pearson" %in% metrics) resultsList$pearson <- cor(x, y, method="pearson")
			if ("rmse" %in% metrics) resultsList$rmse <- sqrt( sum( (x - y)^2 ) / length(x) )
			if ( "meanratio" %in% metrics || subsetStats ) {
				ratio <- y / x
				if ("meanratio"%in% metrics) resultsList$meanratio <- mean(ratio)
				if (subsetStats) resultsList$ratio <- ratio
			}
			# Generate scatter plots
			if (plot && length(x) > 0) {
				# Build output filename prefix for plots
				plotPrefix <- paste(plotPrefix, sampleName, subset, sep=".")
				# Plot
				plotScatter(prefix=plotPrefix, x=x, y=y, xlab=paste(sampleName, "replicate 1", "(log2 expression)"), ylab=paste(sampleName, "replicate 2", "(log2 expression)"), format=format, min=min, max=max)
			}
			# Return results
			return(resultsList)
		})
	}, names(samples), samples, replicates, SIMPLIFY=FALSE)
	# Calculate KS tests if requested
	if (subsetStats) {
		# Extract list of ratios only
		ratios <- lapply(results, lapply, "[[", "ratio")
		# Remove 'ratio' element from results list
		results <- lapply(results, lapply, function(subset) {
			subset$ratio <- NULL
			subset
		})
		# Calculate KS tests & subset results
		ks <- lapply(ratios, calculateKS, statsRefs)
		statsP <- t(sapply(ks, lapply, "[[", "P"))
		statsD <- t(sapply(ks, lapply, "[[", "D"))
		# Build output filenames for tables
		outFileP <- paste(prefix, "ks_P", "tab", sep=".")
		outFileD <- paste(prefix, "ks_D", "tab", sep=".")
		# Write tables
		write.table(statsP, outFileP, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
		write.table(statsD, outFileD, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
	}
	# Get output tables for each metric and write to file
 	if (! metrics == "") {
		invisible(lapply(metrics, function(metric) {
			# Build data table
			tab <- matrix(byrow=TRUE, nrow=length(results), ncol=length(subsets), unlist(sapply(results, sapply, "[", metric)), dimnames=list(names(results), names(subsets) ) )
			# Build output filename for table
			outFile <- paste(prefix, metric, "tab", sep=".")
			# Write table
			write.table(tab, outFile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
		}))
	}
}
#---> CALCULATE METRICS FOR REFERENCE COMPARISONS <---#
calculateReferenceComparisons <- function(samples, references, outDir, outDirPlot, experiment, featureType, subsets, subsetStats, statsRefs, metrics, pseudo, trueFalse, plot, format, min, max) {
	# Iterate over references
	invisible(mapply(function(refName, ref) {
		# Build output filename prefix
		prefix <- file.path(outDir, paste(paste("reference_comparison", refName, sep="_"), experiment, featureType, sep="."))
	        plotPrefix <- file.path(outDirPlot, paste(paste("reference_comparison", refName, sep="_"), experiment, featureType, sep="."))
		# Iterate over samples
		results <- mapply(function(sampleName, sample) {
			# Iterate over subsets
			lapply(names(subsets), function(subset) {
				# Extract sample and reference values
				m <- merge(sample[sample$id %in% subsets[[subset]], ], ref[ref$id %in% subsets[[subset]], ], by=1)
				x <- m[ , 3]	# Reference on the x axis
				y <- m[ , 2]	# Estimates on the y axis
				# Initialize container for results
				resultsList <- list()
				# Calculate metrics
				if ("spearman" %in% metrics) resultsList$spearman <- cor(x, y, method="spearman")
				x[x < pseudo] <- pseudo
				y[y < pseudo] <- pseudo
				x <- log2(x)
				y <- log2(y)
				if ("pearson" %in% metrics) resultsList$pearson <- cor(x, y, method="pearson")
				if ("rmse" %in% metrics) resultsList$rmse <- sqrt(sum((x-y)^2)/length(x))
				if ( "meanratio" %in% metrics || subsetStats ) {
					ratio <- y / x
					if ("meanratio"%in% metrics) resultsList$meanratio <- mean(ratio)
					if (subsetStats) resultsList$ratio <- ratio
				}
				# Generate scatter plots
				if (plot && length(x) > 0) {
					# Build output filename prefix for plots
					plotPrefix <- paste(plotPrefix, sampleName, subset, sep=".")
					# Plot
					plotScatter(prefix=plotPrefix, x=x, y=y, xlab=paste(refName, "(log2 expression)"), ylab=paste(sampleName, "(log2 expression)"), format=format, min=min, max=max)
				}
				# Do type I- / type II-error analysis
				if (trueFalse) {
					expressed <- sum(m[,3] != 0)
					not_expressed <- sum(m[,3] == 0)
					resultsList$true_pos_rate <- sum(m[,3] != 0 & m[,2] != 0) / expressed
					resultsList$true_neg_rate <- sum(m[,3] == 0 & m[,2] == 0) / not_expressed
					resultsList$false_pos_rate <- sum(m[,3] == 0 & m[,2] != 0) / not_expressed
					resultsList$false_neg_rate <- sum(m[,3] != 0 & m[,2] == 0) / expressed
				}
				# Return results
				return(resultsList)
			})
		}, names(samples), samples, SIMPLIFY=FALSE)
		# Calculate KS tests if requested
		if (subsetStats) {
			# Extract list of ratios only
			ratios <- lapply(results, lapply, "[[", "ratio")
			# Remove 'ratio' element from results list
			results <- lapply(results, lapply, function(subset) {
				subset$ratio <- NULL
				subset
			})
			# Calculate KS tests & subset results
			ks <- lapply(ratios, calculateKS, statsRefs)
			statsP <- t(sapply(ks, lapply, "[[", "P"))
			statsD <- t(sapply(ks, lapply, "[[", "D"))
			# Build output filenames for tables
			outFileP <- paste(prefix, "ks_P", "tab", sep=".")
			outFileD <- paste(prefix, "ks_D", "tab", sep=".")
			# Write tables
			write.table(statsP, outFileP, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			write.table(statsD, outFileD, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
		}
		# Get output tables for remaining metrics and write to file
		if (trueFalse) metrics <- c(metrics, "true_pos_rate", "true_neg_rate", "false_pos_rate", "false_neg_rate")
	        if (! metrics == "") {
		 	invisible(lapply(metrics, function(metric) {
				# Build data table
				tab <- matrix(byrow=TRUE, nrow=length(results), ncol=length(subsets), unlist(sapply(results, sapply, "[", metric)), dimnames=list(names(results), names(subsets) ) )
				# Build output filename for table
				outFile <- paste(prefix, metric, "tab", sep=".")
				# Write table
				write.table(tab, outFile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
				return(tab)
			}))
		}
	}, names(references), references))
}
#---> CALCULATE KOLMOGOROV-SMIRNOV STATISTICS <---#
calculateKS <- function(sample, refElement) {
	ref <- sample[[refElement]]
	lapply(sample, function(subset) {
		results <- list()
		test <- suppressWarnings(ks.test(subset, ref))
		results$P <- test$p.value
		results$D <- test$statistic
		return(results)
	})
}
#---> GENERATE SCATTER PLOTS <---#
plotScatter <- function(prefix, x, y, xlab, ylab, format, min, max) {
	if (is.null(min)) min <- floor(min(x))
	if (is.null(max)) max <- ceiling(max(y) * 1.1)
	if ("pdf" %in% format) {
		pdf(paste0(prefix, ".pdf"))
		par(mar=c(5.1,5.1,1.1,1.1))
		heatscatter(x, y, xlim=c(min,max), ylim=c(min,max), xlab=xlab, ylab=ylab, main="", cor=FALSE, ncol=256, cexplot=0.25, cex.axis=1.25, cex.lab=1.25)
		graphics.off()
	}
	if ("svg" %in% format) {
		svg(paste0(prefix, ".svg"))
		par(mar=c(5.1,5.1,1.1,1.1))
		heatscatter(x, y, xlim=c(min,max), ylim=c(min,max), xlab=xlab, ylab=ylab, main="", cor=FALSE, ncol=256, cexplot=0.25, cex.axis=1.25, cex.lab=1.25)
		graphics.off()
	}
	if ("png" %in% format) {
		png(paste0(prefix, ".png"), pointsize=15)
		par(mar=c(5.1,5.1,1.1,1.1))
		heatscatter(x, y, xlim=c(min,max), ylim=c(min,max), xlab=xlab, ylab=ylab, main="", cor=FALSE, ncol=256, cexplot=0.25, cex.axis=1.25, cex.lab=1.25)
		graphics.off()
	}
	if ("png_data_only" %in% format) {
		png(paste0(prefix, ".data_only.png"))
		par(mar=c(0,0,0,0))
		heatscatter(x, y, xlim=c(min,max), ylim=c(min,max), xlab="", ylab="", main="", cor=FALSE, ncol=256, cexplot=0.25, axes=FALSE)
		graphics.off()
	}
}
#================#
#  // FUNCTIONS  #
#================#


#===========#
#  MAIN //  #
#===========#


#---> START MESSAGE <---#
if (verbose) cat("Starting '", script, "'...\n", sep="")


#---> LOAD PACKAGES <---#
if (plot) {

	# Print status message
	if (verbose) cat("Loading packages...\n", sep="")
	
	# Load packages
	if ( suppressWarnings(suppressPackageStartupMessages(require("LSD"))) == FALSE ) { stop("Package 'LSD' required!\nExecution aborted.") }
	if ( suppressWarnings(suppressPackageStartupMessages(require("KernSmooth"))) == FALSE ) { stop("Package 'KernSmooth' required!\nExecution aborted.") }

}

#---> SELECT DATA <---#

# Print status message
if (verbose) cat("Selecting data...\n", sep="")

# Select estimate files
sampleFiles <- dir(path=inDir, pattern=glob2rx(globSample), full.names=TRUE)

# Select replicate files...
if (doReps) {

	sampleFiles <- sampleFiles[basename(sampleFiles) %in% dir(path=repDir)]
	replicateFiles <- file.path(repDir, basename(sampleFiles))

# ...OR reference files
} else referenceFiles <- dir(path=refDir, pattern=glob2rx(globRef), full.names=TRUE)


#---> LOAD DATA <---#

# Print status message
if (verbose) cat("Loading data...\n", sep="")

# Load subsets
load(subsetsFile)
subsetsType <- subsets[[type]]

# Load estimate data
sampleList <- lapply(sampleFiles, read.delim, header=FALSE, stringsAsFactors=FALSE, col.names=c("id", "value"), colClasses=c("character", "numeric"))
names(sampleList) <- unlist(lapply(strsplit(basename(sampleFiles), "\\."), "[", 1))

# Load replicate data...
if (doReps) {

	replicateList <- lapply(replicateFiles, read.delim, header=FALSE, stringsAsFactors=FALSE, col.names=c("id", "value"), colClasses=c("character", "numeric"))
	names(replicateList) <- unlist(lapply(strsplit(basename(sampleFiles), "\\."), "[", 1))

# ...OR reference data
} else {

	referenceList <- lapply(referenceFiles, read.delim, header=FALSE, stringsAsFactors=FALSE, col.names=c("id", "value"), colClasses=c("character", "numeric"))
	names(referenceList) <- unlist(lapply(strsplit(basename(referenceFiles), "\\."), "[", 1))

}


#---> CALCULATE METRICS / STATS & GENERATE SCATTER PLOTS <---#

# Print status message
if (verbose) cat("Calculating metrics...\n", sep="")

# Calculate agreement with replicate...
if (doReps) {
	suppressWarnings(calculateReplicateAgreement(samples=sampleList, replicates=replicateList, outDir=outDir, outDirPlot=outDirPlot, experiment=exp, featureType=type, subsets=subsetsType, subsetStats=subsetStats, statsRefs=statsRefs, metrics=metrics, pseudo=pseudo, plot=plot, format=plotFormat, min=plotMin, max=plotMax))
# ...OR with references
} else {
	suppressWarnings(calculateReferenceComparisons(samples=sampleList, references=referenceList, outDir=outDir, outDirPlot=outDirPlot, experiment=exp, featureType=type, subsets=subsetsType, subsetStats=subsetStats, statsRefs=statsRefs, metrics=metrics, pseudo=pseudo, trueFalse=trueFalse, plot=plot, format=plotFormat, min=plotMin, max=plotMax))
}


#---> END MESSAGE <---#
if (verbose) cat("Done.\n")

#===========#
#  // MAIN  #
#===========#
