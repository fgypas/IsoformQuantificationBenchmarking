#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Apr 30, 2014
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
description <- "Converts the exon entries of a GTF file to a BED12+3 file with one line per transcript.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.1.1 (27-MAY-2015)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--gtf"), action="store", type="character", default="", help="GTF input filename (required; compression via 'gzip' is accepted).", metavar="file"),
		make_option(c("-o", "--bed12"), action="store", type="character", default="", help="BED12 output filename (required).", metavar="file"),
		make_option(c("-g", "--include-gene-id"), action="store_true", default=FALSE, help="Include the gene ID in the name field of the output file. Gene ID and transcript ID will be separated by the delimite specified via --id-delimiter."),
		make_option(c("-d", "--id-delimiter"), action="store", default="$", help="Delimiter to separate gene and transcript ID when --include-gene-id is specified (default: '$'). Ignored if --id-delimiter is not specified.", metavar="delimiter"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die."),
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages to STDOUT.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf <PATH> --bed12 <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if ( opt$gtf == "" || opt$bed12	== "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())	
	stop(print_help(opt_parser))
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
if ( opt$verbose ) cat("Loading packages...\n", sep="")
# Load packages
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$gtf), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> COMPILE LIST OF GENES <---#
# Print status message
if ( opt$verbose ) cat("Compiling transcript list...\n")
# Subset EXONS (discards several other categories, e.g. CDS, start_codon etc.)
gr <- gr[values(gr)[["type"]] == "exon"]
# Merge gene and transcript ID if requested (--include-gene-id option)
if (opt$`include-gene-id`) gr$transcript_id <- paste(gr$gene_id, gr$transcript_id, sep=opt$`id-delimiter`)
# Split exons GRanges into GRangesList by 'transcript_id' 
grl <- split(gr, gr$transcript_id)

#---> EXPORT BED12 <---#
# Print status message
if ( opt$verbose ) cat("Writing output...\n")
# Write output file / adds some track descriptor line to the top that is undesired!
export(object=grl, con=opt$bed12, format="bed15")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
