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
description <- "For each feature or group of features, computes a read count based on a BAM read file.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
created <- "Created: 25-NOV-2014"
version <- "Version: 1.1.1 (27-MAY-2015)"
requirements <- "Requires: optparse, rtracklayer, GenomicAlignments, RSamtools"
msg <- paste(description, author, created, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
                make_option("--bam", action="store", type="character", default="", help="BAM read alignments file (required). A BAI index file must be present, either (1) with the same path and filename plus a '.bai' extension or (2) supplied via '--bai'.", metavar="file"),
                make_option("--bai", action="store", type="character", default="", help="BAM index file (default: argument to '--bam' + '.bai'", metavar="file"),
                make_option("--gtf", action="store", type="character", default="", help="GTF feature file (required).", metavar="file"),
		make_option("--feature-names", action="store", type="character", default="ID", help="Use the specified attribute as the feature name in the output (default: 'ID'). If the selected attribute is not unique, features are grouped and summarized counts are produced.", metavar="attribute"),
                make_option("--counts", action="store", type="character", default="", help="Output filename (required). Table of raw and normalized counts for each feature.", metavar="file"),
                make_option("--normalize-counts", action="store_true", default=FALSE, help="Add columns for counts per million reads (CPM), feature lengths and reads per kilobase per million reads (RPKM)."),
                make_option("--write-rpkm-only", action="store_true", default=FALSE, help="Write a table of only the feature ID and the corresponding RPKM count. Requires '--normalize-counts'."),

                make_option("--header", action="store_true", default=FALSE, help="Add column names."),
                make_option("--mode", action="store", type="character", default="IntersectionStrict", help="Counting mode (one of 'IntersectionStrict', 'Union' or 'IntersectionNotEmpty'; default: 'IntersectionStrict'); refer to the Bioconductor GenomicAlignments manual for details.", metavar="string"),
                make_option("--ignore-strand", action="store_true", default=FALSE, help="Do not consider the strand when counting read/feature overlaps."),
                make_option("--inter-feature", action="store_true", default=FALSE, help="Do not count reads intersecting multiple features."),
                make_option("--paired-end", action="store_false", default=TRUE, help="Specify when counting paired-end reads (default: single-end)."),
                make_option("--include-singletons", action="store_true", default=FALSE, help="Specify if singletons should be included when counting paired-end reads."),
                make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die!"),
                make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die!"),
                make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages to STDOUT.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --bam <PATH> --gtf <PATH> --counts <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if      ( opt$bam == "" || opt$gtf == "" || opt$`feature-names` == "" || opt$counts == "" ) {
        write("[ERROR] Required argument(s) missing!\n\n", stderr())
        stop(print_help(opt_parser))
}
## Set default BAM index filename if not supplied by user
if ( opt$bai == "" ) {
        opt$bai <- paste(opt$bam, "bai", sep=".")
}
## Die if illegal counting mode is specified
if ( ! ( opt$mode == "Union" || opt$mode == "IntersectionStrict" || opt$mode == "IntersectionNotEmpty" ) ) {
        write("[ERROR] Illegal argument to '--mode' option! Please use one of 'Union', 'IntersectionStrict' or 'IntersectionNotEmpty'.\n\n", stderr())
        stop(print_help(opt_parser))
}
## Set '--include-singletons' to FALSE if reads are not paired-ended
if ( opt$`include-singletons` && ! opt$`paired-end`) {
        write("[WARNING] Option '--include-singletons' does not apply when reads are single-ended. Option ignored!\n", stderr())
        opt$`include-singletons` <- FALSE
}
## Set '--write-rpkm-only' to FALSE if '--normalize-counts' is not specified
if ( opt$`write-rpkm-only` && ! opt$`normalize-counts`) {
        write("[WARNING] Option '--write-rpkm-only' does not apply when '--normalize-counts' is not specified. Option ignored!\n", stderr())
        opt$`write-rpkm-only` <- FALSE
}
#===================#
#    OPTIONS END    #
#===================#

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
if ( suppressWarnings(suppressPackageStartupMessages(require("GenomicAlignments"))) == FALSE ) { stop("Package 'GenomicAlignments' required!\nExecution aborted.") }
if ( suppressWarnings(suppressPackageStartupMessages(require("Rsamtools"))) == FALSE ) { stop("Package 'Rsamtools' required!\nExecution aborted.") }

#---> IMPORT GTF <---#
# Print status message
if ( opt$verbose ) cat("Importing features from GTF file '", opt$gtf, "'...\n", sep="")
# Use rtracklayer::import method to import file as GRanges object 
gr <- import(con=opt$gtf, format="gtf", asRangedData=FALSE)

#---> GROUP FEATURES <---#
# Print status message
if ( opt$verbose ) cat("Grouping features by attribute '", opt$`feature-names`, "'...\n", sep="")
# Use gene names for range names
grl <- split(gr, values(gr)[[opt$`feature-names`]])

#---> ACCESS BAM <---#
# Print status message
if ( opt$verbose ) cat("Accessing BAM file '", opt$bam, "'...\n", sep="")
# Use Rsamtools::BamFile and Rsamtools::BamFileList methods to access the file 
bfl <- BamFileList(BamFile(opt$bam, index=opt$bai))

#---> COUNTING OVERLAPS <---#
# Print status message
if ( opt$verbose ) cat("Counting overlaps...\n", sep="")
# Use Rsamtools::BamFile and Rsamtools::BamFileList methods to access the file 
se <- summarizeOverlaps(features=grl, reads=bfl, mode=opt$mode, ignore.strand=opt$`ignore-strand`, inter.feature=opt$`inter-feature`, singleEnd=opt$`paired-end`, fragments=opt$`include-singletons`)
# Build count table
ids <- as.character(rowData(se)@partitioning@NAMES)
raw_counts <- as.integer(assays(se)$counts)
df <- data.frame(id=ids, raw_count=raw_counts)

#---> NORMALIZE READ COUNTS <---#
if ( opt$`normalize-counts` ) {
	# Print status message
	if ( opt$verbose ) cat("Normalizing raw counts...\n", sep="")
	# Extract/calculate additional data
	cpms <- raw_counts / sum(raw_counts) * 1000000
	tmp <- data.frame(id=values(gr)[[opt$`feature-names`]], length=width(gr))
	lengths <- aggregate(length ~ id, tmp, sum)[,2]
	rpkms <- cpms / lengths * 1000
	# Append columns to count table
	df <- cbind(df, cpm=cpms, length=lengths, rpkm=rpkms)
	# Remove all data columns except RPKM values if '--write-rpkm-only' is set
	if ( opt$`write-rpkm-only` ) {
		df <- data.frame(id=df$id, rpkm=df$rpkm)
	}
}

#---> WRITE COUNT TABLE <---#
# Print status message
if ( opt$verbose ) cat("Writing output to table '", opt$counts, "'...\n", sep="")
# Write to output file
write.table(df, file=opt$counts, quote=FALSE, sep="\t", row.names=FALSE, col.names=opt$header)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n\n")

#---> PRINT SESSION INFO <---#
if ( opt$verbose ) {
        cat("++++++++++++\nSESSION INFO\n++++++++++++\n\n")
        writeLines(capture.output(sessionInfo()))
}
#================#
#    MAIN END    #
#================#
