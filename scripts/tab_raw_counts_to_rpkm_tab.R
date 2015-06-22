#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 4, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
#==================#
#    HEADER END    #
#==================#

#===================#
#   OPTIONS START   #
#===================#
#---> LOAD OPTION PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Generates a table of RPKM-normalized read counts from a table of raw read counts and either a table of feature lengths or a GTF gene annotation file.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 04-NOV-2014"
version <- "Version: 1.1.1 (27-MAY-2015)"
requirements <- "Requires: optparse, rtracklayer (only if '--gtf' is specified)"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--raw"), action="store", type="character", default="", help="Raw feature counts in TAB format (required).", metavar="file"),
		make_option(c("-l", "--lengths"), action="store", type="character", default="", help="Gene lengths file in TAB format. Exactly one of '--lengths' and '--gtf' is required.", metavar="file"),
		make_option(c("-g", "--gtf"), action="store", type="character", default="", help="'Union exon' annotations in GTF format (i.e. all exons of a given gene are merged). The sum of the lengths of all 'union exons' of a gene is used for the normalization. Exactly one of '--lengths' and '--gtf' is required.", metavar="file"),
		make_option(c("-o", "--rpkm"), action="store", type="character", default="", help="Normalized (RPKM) counts in TAB format (output filename; required).", metavar="file"),
		make_option(c("-c", "--comprehensive"), action="store_true", default=FALSE, help="Print comprehensive output table including CPM (counts per million reads) and gene/feature lengths."),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die."),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die."),
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --raw <PATH> [--lengths <PATH> --gtf <PATH>] --rpkm <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$raw == "" || (opt$lengths	== "" && opt$gtf == "") || opt$rpkm == "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
	stop(print_help(opt_parser))
}
## Die if mutually exclusive options are specified...
if      ( opt$lengths != "" && opt$gtf != "") {  
        write("[ERROR] Options '--lengths' and '--gtf' are mutually exclusive!\n\n", stderr())
        stop(print_help(opt_parser))
}
#===================#
#    OPTIONS END    #
#===================#

#================#
#   MAIN START   #
#================#
#---> START MESSAGE <---#
if ( opt$verbose ) cat("Starting '", script, "'...\n\n", sep="")

#---> READ TABLE OF RAW COUNTS <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$raw), "'...\n", sep="")
# Read table
raw <- read.table(opt$raw, sep="\t", header=FALSE, col.names=c("id", "raw_count"))

#---> READING/DERIVING GENE/FEATURE LENGHTS <---#
# Print status message
if ( opt$verbose ) cat("Obtain feature/gene lengths for normalization...\n", sep="")
# Either read lengths directly or derive from annotation file
if (opt$lengths != "") {
	#---> READ TABLE OF GENE LENGTHS <---#
	# Print status message
	if ( opt$verbose ) cat("Reading feature lengths from file '", basename(opt$lengths), "'...\n", sep="")
	# Read table
	len <- read.table(opt$lengths, sep="\t", header=FALSE, col.names=c("id", "length"))
} else {
	#---> LOAD PACKAGES <---#
	# Print status message
	if ( opt$verbose ) cat("Loading required packages...\n")
	# Load packages
	if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }
	
	#---> IMPORT GTF <---#
	# Print status message
	if ( opt$verbose ) cat("Reading gene annotation file '", basename(opt$gtf), "'...\n", sep="")
	# Use rtracklayer::import method to import GTF file to GRanges object 
	gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

	#---> SUM UNION EXON LENGTHS <---#
	# Print status message
	if ( opt$verbose ) cat("Deriving 'union exon' gene lengths...\n", sep="")
        # Build a temporary dataframe of gene IDs and 'union exon' lengths
        tmp <- data.frame(id=gr$ID, length=width(gr))
        # Aggregate lengths of individual 'union exons' per gene
        len <- aggregate(length ~ id, tmp, sum)
}

#---> NORMALIZE BY LIBRARY SIZE <---#
# Print status message
if ( opt$verbose ) cat("Compute counts per million reads (CPM)...\n", sep="")
# Compute CPM
cpm <- cbind(raw, cpm=raw$raw_count / sum(raw$raw_count) * 1000000)

#---> NORMALIZE BY GENE LENGTHS <---#
# Print status message
if ( opt$verbose ) cat("Compute reads per kilobase per million reads (RPKM)...\n", sep="")
# Merge CPM table with gene lengths
cpm <- merge(cpm, len, by=1)
# Warn if no size was found for every gene with raw counts
if (nrow(cpm) < nrow(raw)) warning("[WARNING] Gene lengths were not supplied for all genes in the raw count table. Some genes may be missing from the output.")
# Compute RPKM
rpkm <- cbind(cpm, rpkm=cpm$cpm / cpm$length * 1000)

#---> SUBSET TABLE <---#
if ( ! opt$comprehensive ) {
	# Print status message
	if ( opt$verbose ) cat("Subsetting output table...\n", sep="")
	# Subset table
	rpkm <- rpkm[ ,c("id", "rpkm")]
}

#---> WRITE OUTPUT <---#
# Print status message
if ( opt$verbose ) cat("Writing output to TAB file '", opt$rpkm, "'...\n", sep="")
# Write output file
write.table(rpkm, file=opt$rpkm, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
