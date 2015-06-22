#!/usr/bin/Rscript

#=============#
#  HEADER //  #
#=============#
### Created: May 28, 2015
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
description <- "Adds and/or filters features of a normalized expression estimate count table. Gene and transcript expression estimates can be provided separately or individually. If only the latter are provided, the former are determined by aggregation of transcript estimates. Estimates on the level of poly(A) sites are generated likewise, if transcript estimates are available.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 28-MAY-2015"
version <- "Version: 1.0.1 (29-MAY-2015)"
requirements <- "Requires: optparse"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
	make_option("--method", action="store", type="character", default=NULL, help="Method name (required).", metavar="string"),
	make_option("--experiment", action="store", type="character", default=NULL, help="Experiment name (required).", metavar="string"),
	make_option("--output-directory", action="store", type="character", default=NULL, help="Output directory (required).", metavar="path"),
	make_option("--suffix", action="store", type="character", default="estimates", help="Suffix attached to filename (default: 'estimates').", metavar="string"),
	make_option("--no-subdirs", action="store_true", default=FALSE, help="Specify when output files shall be written to a single directory in the form 'output-directory/experiment.feature-type.method.suffix', where 'feature-type' is one of 'transcripts', 'processing_sites' or 'genes' and the other values are derived from the arguments to the options with the same names. By default, output files are written to 'output-directory/experiment/feature-type/method.suffix' (subdirectories are created if not existing)."),
	make_option("--transcript-estimates", action="store", type="character", default=NULL, help="Transcript isoform expression estimates. Exactly one of '--transcript-estimates', '--processing-site-estimates' and '--gene-estimates' is required. Required if '--derive-processing-site-estimates' or '--derive-gene-from-transcript-estimates' is specified.", metavar="path"),
	make_option("--processing-site-estimates", action="store", type="character", default=NULL, help="3'-end/polyA processing site expression estimates. Exactly one of '--transcript-estimates', '--processing-site-estimates' and '--gene-estimates' is required.", metavar="path"),
	make_option("--gene-estimates", action="store", type="character", default=NULL, help="Gene expression estimates. Exactly one of '--transcript-estimates', '--processing-site-estimates' and '--gene-estimates' is required.", metavar="path"),
	make_option("--feature-subsets-R", action="store", type="character", default=NULL, help="R object, output of 'compile_feature_subsets.[hsa|mmu|sim].R' scripts. The feature IDs contained in list element 'all' are used for filtering of the desired feature type. Exactly one of '--feature-subsets-R' and '--id-filter' is required.", metavar="path"),
	make_option("--id-filter", action="store", type="character", default=NULL, help="Text file containing feature IDs to be used for the filtering of the desired feature type. Should contain transcript, processing/polyA site or gene IDs, respectively, when '--transcript-estimates', '--processing-site-estimates' or '--gene-estimates' is specified.  Exactly one of '--feature-subsets-R' and '--id-filter' is required. Options '--id-filter' is not compatible with either of '--derive-processing-site-estimates', '--derive-gene-from-transcript-estimates' and '--derive-gene-from-processing-site-estimates'.", metavar="path"),
	make_option("--derive-processing-site-estimates", action="store_true", default=FALSE, help="In addition to the specified estimates, derive 3'-end/polyA processing site expression from transcript/isoform expression estimates by summing the expression estimates of all transcripts ending in a given processing site. Requires that a transcript ID to processing site ID lookup table is supplied as an argument to option '--transcript-to-processing-site-id'. Options '--derive-processing-site-estimates' and '--id-filter' are mutually exclusive."),
	make_option("--derive-gene-from-transcript-estimates", action="store_true", default=FALSE, help="In addition to the specified estimates, derive gene expression from transcript/isoform expression estimates by summing the expression estimates of all transcripts annotated for a given gene. Requires that a transcript ID to gene ID lookup table is supplied as an argument to option '--transcript-to-gene-id'. Options '--derive-gene-from-transcript-estimates', '--derive-gene-from-processing-site-estimates' and '--id-filter' are mutually exclusive."),
	make_option("--derive-gene-from-processing-site-estimates", action="store_true", default=FALSE, help="In addition to the specified estimates, derive gene expression from 3'-end/polyA processing site expression estimates by summing the expression estimates of all processing sites annotated for a given gene. Requires that a processing site ID to gene ID lookup table is supplied as an argument to option '--processing-site-to-gene-id'. Options '--derive-gene-from-processing-site-estimates', '--derive-gene-from-transcript-estimates' and '--id-filter' are mutually exclusive."),
	make_option("--transcript-to-processing-site-id", action="store", type="character", default=NULL, help="Tab-separated lookup table of the form 'Transcript ID -> 3'-end/polyA processing site ID'. Required if '--derive-processing-site-estimates' is specified.", metavar="path"),
	make_option("--transcript-to-gene-id", action="store", type="character", default=NULL, help="Tab-separated lookup table of the form 'Transcript ID -> gene ID'. Required if '--derive-gene-from-transcript-estimates' is specified.", metavar="path"),
	make_option("--processing-site-to-gene-id", action="store", type="character", default=NULL, help="Tab-separated lookup table of the form '3'-end/polyA processing site ID -> gene ID'. Required if '--derive-gene-from-processing-site-estimates' is specified.", metavar="path"),
	make_option("--help", action="store_true", default=FALSE, help="Show this information and die."),
	make_option("--usage", action="store_true", default=FALSE, dest="help", help="Show this information and die."),
	make_option("--verbose", action="store_true", default=FALSE, help="Print log messages.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --method <STRING> --experiment <STRING> --output-directory <PATH> [--transcript-estimates <PATH> --processing-site-estimates <PATH> --gene-estimates <PATH>]  [--feature-subsets-R <PATH> --id-filter <PATH>]\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Re-assign command-line arguments
meth <- opt$`method`
exp <- opt$`experiment`
outDir <- opt$`output-directory`
suffix <- opt$`suffix`
noSubDirs <- opt$`no-subdirs`
estTrx <- opt$`transcript-estimates`
estPas <- opt$`processing-site-estimates`
estGen <- opt$`gene-estimates`
subsetsR <- opt$`feature-subsets-R`
filter <- opt$`id-filter`
getPas <- opt$`derive-processing-site-estimates`
getGenFromTrx <- opt$`derive-gene-from-transcript-estimates`
getGenFromPas <- opt$`derive-gene-from-processing-site-estimates`
trxToPas <- opt$`transcript-to-processing-site-id`
trxToGen <- opt$`transcript-to-gene-id`
pasToGen <- opt$`processing-site-to-gene-id`
verbose <- opt$`verbose`

## Die if required arguments are missing
if ( is.null(meth) || is.null(exp) || is.null(outDir) ) {
        write("[ERROR] Required argument(s) missing!\n\n", stderr())
        stop(print_help(opt_parser))
}
## Die if none or more than one option of a set of mutually exclusive, yet required options are specified
if ( ! sum( c(! is.null(estTrx), ! is.null(estPas), ! is.null(estGen) ) ) == 1 ) {
        write("[ERROR] Mutually exclusive, yet required options are specified or missing! Specify exactly one of '--transcript-estimates', '--processing-site-estimates' and '--gene-estimates'.\n\n", stderr())
        stop(print_help(opt_parser))
}
if ( ! sum( c(! is.null(subsetsR), ! is.null(filter) ) ) == 1 ) {
        write("[ERROR] Mutually exclusive, yet required options are specified or missing! Specify exactly one of '--feature-subsets-R' and '--id-filter'.\n\n", stderr())
        stop(print_help(opt_parser))
}
## Die if mutually exclusive options are specified
if ( sum( c(getGenFromTrx, getGenFromPas) ) > 1 ) {
        write("[ERROR] Mutually exclusive options specified! Specify at most one of '--derive-gene-from-transcript-estimates' and '--derive-gene-from-processing-site-estimates'.\n\n", stderr())
        stop(print_help(opt_parser))
}
if ( ! is.null(filter) && sum( c(getPasFromTrx, getGenFromTrx, getGenFromPas) ) > 0 ) {
        write("[ERROR] Mutually exclusive options specified! Option '--id-filter' is not compatible with either of '--derive-processing-site-estimates', '--derive-gene-from-transcript-estimates' and '--derive-gene-from-processing-site-estimates'.\n\n", stderr())
        stop(print_help(opt_parser))
}
## Die if dependent options are not specified
if ( getPas && ( is.null(estTrx) || is.null(trxToPas) ) ) {
        write("[ERROR] Required argument(s) missing! When specifying option '--derive-processing-site-estimates', options '--transcript-estimates' and '--transcript-to-processing-site-id' need to be specified!\n\n", stderr())
        stop(print_help(opt_parser))
}
if ( getGenFromTrx && ( is.null(estTrx) || is.null(trxToGen) ) ) {
        write("[ERROR] Required argument(s) missing! When specifying option '--derive-gene-from-transcript-estimates', options '--transcript-estimates' and '--transcript-to-gene-id' need to be specified!\n\n", stderr())
        stop(print_help(opt_parser))
}
if ( getGenFromPas && ( is.null(pasToGen) || ( is.null(estPas) && ( is.null(estTrx) || is.null(trxToPas) ) ) ) ) {
        write("[ERROR] Required argument(s) missing! When specifying option '--derive-gene-from-processing-site-estimates', option '--processing-site-to-gene-id' and either '--processing-site-estimates' or both '--transcript-estimates' and '--transcript-to-processing-site-id' need to be specified!\n\n", stderr())
        stop(print_help(opt_parser))
}
#==============#
#  // OPTIONS  #
#==============#


#================#
#  FUNCTIONS //  #
#================#
#---> AGGREGATE EXPRESSION ESTIMATES <---#
aggregateExpression <- function(estimates, lookup_table) {
	df <- merge(estimates[ , 1:2], lookup_table[ , 1:2])
	df <- aggregate(df[ , 2] ~ df[ , 3], df, sum)
	colnames(df) <- c(colnames(lookup_table)[1], colnames(estimates)[2])
	return(df)
}
#---> BUILD OUTPUT FILENAME <---#
buildOutputFilename <- function(dir, pre1, pre2, main, suff, noSubDirTree=FALSE, sep=".") {
	if (noSubDirTree) {
		file <- paste(pre1, pre2, main, suff, sep=sep)
	} else {
		dir <- file.path(dir, pre1, pre2)
		file <- paste(main, suff, sep=sep)
	}
	dir.create(dir, recursive=TRUE, showWarnings=FALSE)
	file <- file.path(dir, file)
	return(file)
}
#================#
#  // FUNCTIONS  #
#================#


#===========#
#  MAIN //  #
#===========#

#---> START MESSAGE <---#
if (verbose) cat("Starting '", script, "'...\n", sep="")


#---> LOAD RESOURCES <---#

# Print status message
if (verbose) cat("Loading resources...\n", sep="")

# Load lookup tables
if ( ! is.null(trxToPas) ) trx2pas <- read.table(trxToPas, sep="\t", header=FALSE, col.names=c("id", "pas_id"), colClasses=c("character", "character"))
if ( ! is.null(trxToGen) ) trx2gen <- read.table(trxToGen, sep="\t", header=FALSE, col.names=c("id", "gen_id"), colClasses=c("character", "character"))
if ( ! is.null(pasToGen) ) pas2gen <- read.table(pasToGen, sep="\t", header=FALSE, col.names=c("id", "gen_id"), colClasses=c("character", "character"))

# Load feature subsets or ID filter
if ( ! is.null(subsetsR) ) load(subsetsR)
if ( ! is.null(filter) ) filt <- scan("tmpfile", what="character", sep="\n")


#---> LOAD ESTIMATES <---#

# Initialize output container and load estimates
outputs <- list()

# If transcript estimates are provided...
if ( ! is.null(estTrx) ) {

	# Print status message
	if (verbose) cat("Loading estimates from file '", estTrx, "'...\n", sep="")

	# Load estimates
	outputs$trx  <- read.table(estTrx, sep="\t", header=FALSE, col.names=c("id", "value"), colClasses=c("character", "numeric"))

}

# If processing site estimates are provided...
if ( ! is.null(estPas) ) {

	# Print status message
	if (verbose) cat("Loading estimates from file '", estPas, "'...\n", sep="")

	# Load estimates
	outputs$pas  <- read.table(estPas, sep="\t", header=FALSE, col.names=c("id", "value"), colClasses=c("character", "numeric"))

}

# If gene estimates are provided...
if ( ! is.null(estGen) ) {

	# Print status message
	if (verbose) cat("Loading estimates from file '", estGen, "'...\n", sep="")

	# Load estimates
	outputs$gene <- read.table(estGen, sep="\t", header=FALSE, col.names=c("id", "value"), colClasses=c("character", "numeric"))

}


#---> AGGREGATE ESTIMATES <---#

# If estimate aggregation is requested...
if ( getPas || getGenFromTrx || getGenFromPas ) {

	# Print status message
	if (verbose) cat("Aggregating estimates...\n", sep="")

	# Get processing site expression from transcript expression
	if ( getPas ) {
		outputs$pas  <- aggregateExpression(outputs$trx, trx2pas)
	}

	# Get gene expression from transcript expression
	if ( getGenFromTrx ) {
		outputs$gene <- aggregateExpression(outputs$trx, trx2gen)
	}

	# Get gene expression from processing site expression	
	if ( getGenFromPas ) {
		outputs$gene <- aggregateExpression(outputs$pas, pas2gen)
	}

}


#---> PROCESS DATA <---#

# Print status message
if (verbose) cat("Processing data...\n", sep="")

# Iterate over output container
outputs <- mapply(function(estimates, type) {

	# Get reference ID set
	if ( is.null(filter) ) {
		reference <- subsets[[type]]$all
	} else {
		reference <- filt
	}

	# Filter features not in reference ID set
        estimates <- estimates[estimates[ , 1] %in% reference, ]

	# Substitute NAs with zeros
	estimates$value[is.na(estimates$value)] <- 0

	# Add features from reference ID set for which no estimates are present (values added: zeros)
	ids_to_add <- reference[! reference %in% estimates[ , 1]]
	values_to_add <- rep(0, length(ids_to_add))
	estimates <- rbind(estimates, data.frame(id=ids_to_add, value=values_to_add))

	# Return estimates
	return(estimates)

}, outputs, names(outputs), SIMPLIFY=FALSE)


#---> WRITE OUTPUT <---#

# Iterate over output container
dump <- mapply(function(estimates, type) {

	# Build output filename
	outFile <- buildOutputFilename(dir=outDir, pre1=exp, pre2=type, main=meth, suff=suffix, noSubDirTree=noSubDirs)

        # Print status message
        if (verbose) cat("Writing output to file '", outFile, "'...\n", sep="")
	
	# Write estimates table
	write.table(estimates, outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

}, outputs, names(outputs))


#---> END MESSAGE <---#
if (verbose) cat("Done.\n")

#===========#
#  // MAIN  #
#===========#
