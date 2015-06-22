#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: sam_split_body_header.pl
### Created: Aug 28, 2013
### Modified: Aug 04, 2014
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: GetOpt::Long
#==================#
### Description: Splits a SAM file and puts out either its header, its body or both. It is assumed that the header is located in one piece at the top of the input file.
### Output: The body of a SAM file, its header or both
### Usage: perl ./sam_split_body_header.pl for information on required and optional arguments
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;


#---> USAGE <---#
my $usage_info = &usage;

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = '';
my $sam = '';
my $body = '';
my $head = '';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	#-----------------------#
	'sam=s' => \$sam,
	'body=s' => \$body,
	'head=s' => \$head
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$sam || !($head || $body); 
die "[ERROR] Could not find input file '$sam'.\n" unless -e $sam;

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'sam_split_body_header.pl'...\n" unless $quiet;

#---> BODY <---#
&sam_to_header_and_body($sam, $body, $head);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit 0;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current script
### Accepts: n/a
### Returns: String with usage information
### Type: Specialized
'Usage: perl ./sam_split_body_header.pl [OPTIONS] --sam [SAM FILE] (--body [SAM FILE BODY] --head [SAM FILE HEADER])

Description: Extracts the header and/or the body of a SAM file. It is assumed that the header is in accordance with the SAM format specifications and is located at the beginning of the document in one continuous stretch.

--sam	Input file in SAM format [REQUIRED]
--body	Body of input SAM file (at least one of --body or --head is required)
--head	Header of input SAM file (at least one of --body or --head is required)
--help	Show this information and die
--quiet	Shut up!
';
}
#-----------------------#
sub sam_to_header_and_body {
### Function: function
### Accepts: 1. Input file in SAM format; 2. Output file for SAM body (or empty string); 3. Output file SAM header (or empty string)
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($sam, $body, $head) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Splitting SAM file '${sam}...'\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $regex_header = '^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$';
	my $regex_comment = '/^\@CO\t.*/';
	my $last_line;
	
	#---> BODY <---#
	
		## Open file handles SAM, BODY and HEAD
		open SAM, "<", $sam or die "[ERROR] Could not open file '$sam'!\n";
		open BODY, ">", $body or die "[ERROR] Could not open file '$body'!\n" if $body;
		open HEAD, ">", $head or die "[ERROR] Could not open file '$head'!\n" if $head;
		
		#---> Check if lines are header lines, then print to dedicated file if requested; once non-header line is met, exit loop <---# 
		## Traverse through file handle SAM line by line
		while (<SAM>) {
			## If line is a header line, print to filehandle HEAD (if $head is set)
			if ( /$regex_header/ || /$regex_comment/ ) {
				print HEAD if $head;
			}
			## Else exit loop
			else {
				$last_line = $_;
				last;
			}
		}
	
		#---> Print body (i.e. rest of lines) to dedicated file if requested <---#
		## Check whether $body is set
		if ( $body ) {
			## If so, set the correct filehandle position...
			{
				use bytes;
				seek SAM, -length($last_line), 1;			
			}
			## ... and print all the rest of filehandle SAM to filehandle BODY
			while (<SAM>) {
				print BODY;
			}
		}
	
		## Close file handles SAM, BODY and HEAD
		close HEAD if $head;
		close BODY if $body;
		close SAM;
	
	#---> STATUS MESSAGES <---#
	print STDOUT "SAM file split." . "\n" unless $quiet;
	print STDOUT "SAM header written to file '${head}'.\n" if $head;
	print STDOUT "SAM body written to file '${body}'.\n" if $body;
	
	#---> RETURN VALUE <---#
	return ;
}
#=======================#
#    SUBROUTINES END    #
#=======================#
