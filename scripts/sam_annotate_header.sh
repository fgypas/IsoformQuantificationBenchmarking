#!/bin/bash

#==================#
#   HEADER START   #
#==================#
### Name: sam_annotate_header.sh
### Created: Feb 17, 2013
### Modified: Feb 17, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: Bash (not POSIX compliant)
#==================#
#    HEADER END    #
#==================#

# Define variable containing usage information
script=`basename $0`
usage="bash ./$script sam_genome_header [FILE] comment_lines [FILE]"
descr='Description:\nAppends comment lines to the end of a SAM header file. Only valid comment lines of the form "@CO \tab\ COMMENT..." are added, all other lines are ignored.'

## Die if wrong number of arguments
if [ $# -ne 2 ]; then
	echo -e "[ERROR] Wrong number of arguments!\n\nUsage: $usage\n\n${descr}\n" >&2
	exit 1
fi

## Die if input files do not exist
for file in $*; do
	if [ ! -f $file ]; then
		echo -e "[ERROR] File '$file' not found!\n\nUsage: $usage\n\n${descr}\n" >&2
		exit 1
	fi
done

# Only proper comment lines are used
com=`perl -n -e 'print if /^\@CO\t/' < $2`

# Merge header and comments
echo -e "$com" | cat $1 -

# Exit
exit 0