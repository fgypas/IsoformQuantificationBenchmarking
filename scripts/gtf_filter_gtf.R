#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 24, 2014
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD OPTION PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "<Extracts features of specified types from a GTF file and optionally filters them by fields 'transcript_id' and 'transcript_type'.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0.1 (20-MAY-2015)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
                make_option(c("-i", "--gtf-in"), action="store", type="character", default="", help="GTF input filename (required).", metavar="file"),
                make_option(c("-o", "--gtf-out"), action="store", type="character", default="", help="GTF output filename (required).", metavar="file"),
                make_option(c("-f", "--feature-type-filter"), action="store", type="character", default="transcript", help="Names of one ore more feature types. Only entries corresponding to these types will be kept in the output (default: 'transcript'). Multiple entries must be separated by comma and no white space in between.", metavar="file"),
                make_option(c("-n", "--transcript-id-filter"), action="store", type="character", default="", help="File containing transcript IDs. Only entries corresponding to these IDs will be kept in the output (default: all entries are kept). IDs must be supplied in a file with one ID per line or in a column of a tab-delimited file (see option '--column').", metavar="file"),
                make_option(c("-c", "--column"), action="store", type="integer", default=1, help="Index (1-based) of column containing the transcript IDs used for filtering (default: 1). Only relevant if an argument to option '--filter' was supplied.", metavar="int"),
                make_option(c("-t", "--transcript-type-filter"), action="store", type="character", default="", help="Names of one ore more transcript types. Only entries corresponding to these types will be kept in the output (default: all entries are kept). Multiple entries must be separated by comma and no white space in between.", metavar="file"),
                make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die."),
                make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die."),
                make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages to STDOUT.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf-in <PATH> --gtf-out <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if ( opt$`gtf-in` == "" || opt$`gtf-out` == "" ) {
        write("[ERROR] Required argument(s) missing!\n\n", stderr())
        stop(print_help(opt_parser))
}
## Get vector of feature types
if ( ! opt$`feature-type-filter` == "transcript" ) {
        opt$`feature-type-filter` <- unlist(strsplit(opt$`feature-type-filter`, split=","))
}
## Get vector of transcript types
if ( ! opt$`transcript-type-filter` == "" ) {
	opt$`transcript-type-filter` <- unlist(strsplit(opt$`transcript-type-filter`, split=","))
}
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n", sep="")

#---> LOAD PACKAGES <---#
# Print status message
if ( opt$verbose ) cat("Loading required packages...\n")
# Load packages
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading GTF file '", basename(opt$`gtf-in`), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$`gtf-in`, format="gtf", asRangedData=FALSE)

#---> FILTER FEATURES <---#
# Print status message
if ( opt$verbose ) cat("<Extracting selected feature type(s)...\n")
# Subset selected features
gr <- gr[gr$type %in% opt$`feature-type-filter`]
# Test if at least one region is returned
if ( ! is.null(nrow(gr)) ) stop("No entries of specified feature type(s) in input file! Check file.\nExecution aborted.\n")

#---> SUBSET FEATURES BY TRANSCRIPT ID <---#
if ( ! opt$`transcript-id-filter` == "" ) {

	#---> LOAD IDS TO SUBSET <---#
	# Print status message
	if ( opt$verbose ) cat("Reading transcript IDs from file '", basename(opt$`transcript-id-filter`), "'...\n", sep="")
	## Read transcript IDs
	subset <- read.table(opt$`transcript-id-filter`, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	if ( ncol(subset) < opt$column ) stop("The specified column index is out of bounds! Check the column index and make sure the filter file is tab-delimited.\nExecution aborted.") 
	subset <- as.character(subset[ , opt$column, drop=TRUE])
	
	#---> SUBSET FEATURES BY IDS <---#
	# Print status message
	if ( opt$verbose ) cat("Filtering features by supplied transcript IDs...\n")
	# Subset
	gr <- gr[gr$transcript_id %in% subset]
}

#---> SUBSET FEATURES BY TRANSCRIPT TYPE <---#
if ( ! opt$`transcript-type-filter`[[1]] == "" ) {

        #---> SUBSET FEATURES BY IDS <---#
        # Print status message
        if ( opt$verbose ) cat("Filtering features by supplied transcript types...\n")
        # Subset
        gr <- gr[gr$transcript_type %in% opt$`transcript-type-filter`]
}

#---> EXPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Writing features to GTF file '", opt$`gtf-out` , "'...\n", sep="")
# Export GTF file
export(object=gr, con=opt$`gtf-out`, format="gtf")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
