#!/usr/bin/Rscript

#==================#
#   HEADER START //  #
#==================#
### Created: May 21, 2015
### Author: Foivos Gypas
### Company: Zavolan Group, Biozentrum, University of Basel
#==================#
#  //  HEADER END    #
#==================#

#===================#
#   OPTIONS START   #
#===================#
#---> LOAD OPTION PARSER <---#
if ( suppressWarnings(suppressPackageStartupMessages(require("optparse"))) == FALSE ) { stop("Package 'optparse' required!\nExecution aborted.") }

#---> GET SCRIPT NAME <---#
script <- sub("--file=", "", basename(commandArgs(trailingOnly=FALSE)[4]))

#---> DESCRIPTION <---#
description <- "Exponantiates the values in a column of a tab-separated file."
author <- "Author: Foivos Gypas & Alexander Kanitz, Biozentrum, University of Basel"
version <- "Version: 1.0.1 (27-MAY-2015)"
requirements <- "Requires: optparse"
msg <- paste(description, author, version, requirements, sep="\n")

#---> COMMAND-LINE ARGUMENTS <---#
## List of allowed/recognized arguments
option_list <- list(
                make_option(c("-i", "--input"), action="store", type="character", default="", help="Input of MMSEQ estimation (logmu; required).", metavar="file"),
                make_option(c("-o", "--output"), action="store", type="character", default="", help="Output filename (required).", metavar="file"),
                make_option(c("-c", "--column"), action="store", type="integer", default=2, help="Column (1-based) that contains values to exponentiate (default: 2).", metavar="file"),
                make_option(c("-h", "--help"), action="store_true", default=FALSE, help="Show this information and die."),
		make_option(c("-u", "--usage"), action="store_true", default=FALSE, dest="help", help="Show this information and die.")
)

## Parse command-line arguments
opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS] --input <PATH> --output <PATH>\n", option_list = option_list, add_help_option=FALSE, description=msg)
opt <- parse_args(opt_parser)

## Die if any required arguments are missing...
if ( opt$`input` == "" || opt$`output` == "" ) {
        write("[ERROR] Required argument(s) missing!\n\n", stderr())
        stop(print_help(opt_parser))
}


#================#
#   MAIN START   #
#================#

# Load data
df <- read.table(opt$input, header=FALSE, sep="\t")

# Exponentiate values in selected column
df[opt$column] <- exp(df[opt$column])

# Write data
write.table(df, opt$output, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

#================#
#    MAIN END    #
#================#
