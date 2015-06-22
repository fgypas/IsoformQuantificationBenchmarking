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
description <- "Extracts the exon entries of a GTF file and optionally filters them by fields 'transcript_id' and 'transcript_type'.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0.1 (27-MAY-2015)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
                make_option(c("-i", "--gtf-in"), action="store", type="character", default="", help="GTF input filename (required).", metavar="file"),
                make_option(c("-o", "--gtf-out"), action="store", type="character", default="", help="GTF output filename (required).", metavar="file"),
                make_option(c("-n", "--name-filter"), action="store", type="character", default="", help="File containing transcript IDs. Only entries corresponding to these IDs will be kept in the output (default: all entries are kept). IDs must be supplied in a file with one ID per line or in a column of a tab-delimited file (see option '--column').", metavar="file"),
                make_option(c("-c", "--column"), action="store", type="integer", default=1, help="Index (1-based) of column containing the transcript IDs used for filtering (default: 1). Only relevant if an argument to option '--filter' was supplied.", metavar="int"),
                make_option(c("-t", "--type-filter"), action="store", type="character", default="", help="Name of one ore more transcript types. Only entries corresponding to these types will be kept in the output (default: all entries are kept). Multiple entries must be separated by comma and no white space in between.", metavar="file"),
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
## Get vector of transcript types
if ( ! opt$`type-filter` == "" ) {
	opt$`type-filter` <- unlist(strsplit(opt$`type-filter`, split=","))
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

#---> FILTER EXONS <---#
# Print status message
if ( opt$verbose ) cat("Extracting entries of type 'exon'...\n")
# Subset exonS (discard all other categories, e.g. CDS, start_codon etc.)
gr <- gr[gr$type == "exon"]
# Test if at least one region is returned
if ( ! is.null(nrow(gr)) ) stop("No entries of type 'exon' in input file! Check file.\nExecution aborted.\n")

#---> SUBSET EXONS BY TRANSCRIPT ID <---#
if ( ! opt$`name-filter` == "" ) {

	#---> LOAD IDS TO SUBSET <---#
	# Print status message
	if ( opt$verbose ) cat("Reading transcript IDs from file '", basename(opt$`name-filter`), "'...\n", sep="")
	## Read transcript IDs
	subset <- read.table(opt$`name-filter`, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	if ( ncol(subset) < opt$column ) stop("The specified column index is out of bounds! Check the column index and make sure the filter file is tab-delimited.\nExecution aborted.") 
	subset <- as.character(subset[ , opt$column, drop=TRUE])
	
	#---> SUBSET EXONS BY IDS <---#
	# Print status message
	if ( opt$verbose ) cat("Filtering exons by supplied transcript IDs...\n")
	# Subset
	gr <- gr[gr$transcript_id %in% subset]
}

#---> SUBSET EXONS BY TRANSCRIPT TYPE <---#
if ( ! opt$`type-filter`[[1]] == "" ) {

        #---> SUBSET EXONS BY IDS <---#
        # Print status message
        if ( opt$verbose ) cat("Filtering exons by supplied transcript types...\n")
        # Subset
        gr <- gr[gr$transcript_type %in% opt$`type-filter`]
}

#---> EXPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Writing exons to GTF file '", opt$`gtf-out` , "'...\n", sep="")
# Export GTF file
export(object=gr, con=opt$`gtf-out`, format="gtf")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#

