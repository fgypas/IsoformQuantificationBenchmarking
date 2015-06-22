#!/usr/bin/perl

#=============#
#  HEADER //  #
#=============#
## Created: Jul 16, 2014
## Author: Alexander Kanitz
## Company: Zavolan Group, Biozentrum, University of Basel
## Requirements: Getopt::Long
#=============#
#  // HEADER  #
#=============#


#========================#
#  PRAGMAS & MODULES //  #
#========================#
use strict;
use warnings;
use Getopt::Long;


#======================#
#  USAGE & VERSION //  #
#======================#
sub usage {
### Returns usage information for current script in a string
<<USAGE;
Usage: perl $0 <FASTA>

Description: Computes transcript lengths from one or more FASTA files and prints them to STDOUT.

Options:
        --trim-id               From each ID, remove everything following the first whitespace (default: do not trim).
        --usage | --help        Show this screen and exit.
        --version               Show version information and exit.
        --verbose               Print log information to STDERR.

Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Jul 16, 2014.
Version 1.0.1 (May 20, 2015)
USAGE
}
#-----------------------#
sub version {
### Returns version information for current script in a string
<<VERSION;
$0 version 1.0.1 (May 20, 2015)
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Jul 16, 2014.
VERSION
}
#======================#
#  // USAGE & VERSION  #
#======================#


#====================#
#  PARSE OPTIONS //  #
#====================#
#---> OPTIONS VARIABLES <---#
my $trim_id = 0;
my $usage = 0;
my $version = 0;
my $verbose= 0;

#---> PARSE / ASSIGN OPTIONS <---#
my $options_result = GetOptions (
        'trim-id' => \$trim_id,
        'usage|help' => \$usage,
        'version' => \$version,
        'verbose' => \$verbose
);

#---> VERIFY OPTIONS <---#
# Print usage information and exit if option parsing was not successful
die &usage unless $options_result;

# Print usage information and exit if --usage or --help were specified
die &usage if $usage;

# Print version information and exit if --version was specified
die &version if $version;
#====================#
#  // PARSE OPTIONS  #
#====================#


#===========#
#  MAIN //  #
#===========#
#---> STATUS MESSAGE <---#
print STDERR "Starting '$0'...\n" if $verbose;

#---> MAIN VARIABLES <---#
my $seq;

        #---> BODY <---#

        #---> Iterate over lines of input file(s) <---#
	while (<>)
	{

		#---> Remove trailing newline character <---#
		chomp;

		#---> IF line is header line... <---#
		if (s/^>//) {

			#---> Print length of previous sequence unless no previous sequence has been encountered (first entry) <---#
			print length($seq), "\n" unless not defined $seq;

			#---> Reset sequence string <---#
			$seq = undef;
			
			#---> Get ID of new sequence (string between '>' and newline or, if '--trim-id' is set, up until first 'whitespace') <---#
			my $id = $_;
			if ($trim_id) {
				$id = (split /\s/, $id, 2)[0];
			}

			#---> Print ID <---#
                        print $id, "\t";

		}

		#---> ...ELSE append current sequence line (excluding the newline) to existing sequence string
		else { $seq .= $_; }
	}

	#---> Complete last entry: Print length of sequence unless no previous sequences has been encountered (empty/invalid file) <---#
	print length($seq), "\n" unless not defined $seq;

#---> STATUS MESSAGE <---#
print STDERR "Done.\n" if $verbose;

#---> PROGRAM EXIT <---#
exit 0;
#===========#
#  // MAIN  #
#===========#
