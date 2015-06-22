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


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> LOAD OPTION PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Removes ambiguous associations (non-unique keys) from a lookup table of the format KEY (tab) VALUE. All entries containing such keys are discarded.\n"
author <- "Author: Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0.1 (20-MAY-2015)"
requirements <- "Requires: optparse"
msg <- paste(description, author, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
                make_option(c("-i", "--infile"), action="store", type="character", default="", help="Tab-separated file of the format KEY (tab) VALUE.", metavar="file"),
                make_option(c("-o", "--outfile"), action="store", type="character", default="", help="Output file (same format as input).", metavar="file"),
                make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die."),
                make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die."),
                make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print log messages to STDOUT.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --infile <PATH> --outfile <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if ( opt$`infile` == "" || opt$`outfile` == "" ) {
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

#---> READ DATA <---#
# Print status message
if ( opt$verbose ) cat("Reading input file '", opt$infile, "'...\n", sep="")
# Read input file
df <- read.table(opt$infile, sep="\t", header=FALSE, col.names=c("key", "value"), colClasses=c("character", "character"))

#---> FIND AMBIGUOUS ASSOCIATIONS <---#
# Print status message
if ( opt$verbose ) cat("Searching for non-unique keys...\n")
# Find non-unique keys
ambiguous <- table(df$key)
ambiguous <- names(ambiguous)[ambiguous > 1]

#---> DISCARD AMBIGUOUS ASSOCIATIONS <---#
# Print status message
if ( opt$verbose ) cat("Discarding entries with non-unique keys...\n")
# Remove entries with non-unique keys
df <- df[! df$key %in% ambiguous, ]

#---> WRITE TABLE <---#
# Print status message
if ( opt$verbose ) cat("Writing lookup table to file '", opt$outfile , "'...\n", sep="")
# Write table
write.table(df, opt$outfile, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

#---> END MESSAGE <---#
if ( opt$verbose ) cat("Done.\n")
#================#
#    MAIN END    #
#================#
