#!/usr/bin/Rscript

#==================#
#   HEADER START   #
#==================#
### Created: Nov 25, 2014
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
description <- "<From a GTF gene set annotation file, generates a GTF file containing the 'pseudoexons' resulting from the union of all exons of all transcripts of each gene (one row per gene).\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 25-NOV-2014"
version <- "Version: 1.1.1 (27-MAY-2015)"
requirements <- "Requires: optparse, rtracklayer"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
		make_option(c("-i", "--gtf-in"), action="store", type="character", default="", help="GTF input filename (required).", metavar="file"),
		make_option(c("-o", "--gtf-out"), action="store", type="character", default="", help="GTF output filename (required).", metavar="file"),
		make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die!"),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die!"),
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages to STDOUT.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --gtf-in <PATH> --gtf-out <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if 	( opt$`gtf-in`	== ""	||	opt$`gtf-out`	== "" ) { 
	write("[ERROR] Required argument(s) missing!\n\n", stderr())
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

#---> LOAD PACKAGES <---#
# Print status message
if ( opt$verbose ) cat("Loading required packages...\n")
# Load packages
if ( suppressWarnings(suppressPackageStartupMessages(require("rtracklayer"))) == FALSE ) { stop("Package 'rtracklayer' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", basename(opt$`gtf-in`), "'...\n", sep="")
# Use rtracklayer::import method to import GTF file to GRanges object 
gr <- import(con=opt$`gtf-in`, format="gtf", asRangedData=FALSE)

#---> SUBSET EXONS <---#
# Print status message
if ( opt$verbose ) cat("Subsetting exons...\n")
# Subset EXONS (discards all other categories, e.g. CDS, start_codon etc.)
gr <- gr[values(gr)[["type"]] == "exon"]
# Test if at least one region is returned
if ( length(gr) == 0 ) stop("No entries of type 'exon' in input file! Check file.\nExecution halted.\n")

#---> COMPILE LIST OF GENES <---#
# Print status message
if ( opt$verbose ) cat("Grouping exons by genes...\n")
# Split exons GRanges into GRangesList by 'gene_id' 
grl <- split(gr, gr$gene_id)

#---> TAKE UNION OF EXONS <---#
# Print status message
if ( opt$verbose ) cat("Merging exons...\n")
# Make union of exons to generate pseudoexons
grl <- reduce(grl, min.gapwidth=1L)

#---> UNLIST GROUPED PSEUDOEXONS <---#
# Print status message
if ( opt$verbose ) cat("Ungrouping 'pseudoexons'...\n")
# Unlist GRangesList
gr <- unlist(grl)
	
#---> EXPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Writing output to GTF file '", opt$`gtf-out`, "'...\n", sep="")
# Write output file
export(object=gr, con=opt$`gtf-out`, format="gtf")

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
