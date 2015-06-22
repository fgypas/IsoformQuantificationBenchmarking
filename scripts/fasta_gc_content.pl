#!/usr/bin/perl

#=============#
#  HEADER //  #
#=============#
## Created: Sep 22, 2014
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

#========================#
#  // PRAGMAS & MODULES  #
#========================#


#======================#
#  USAGE & VERSION //  #
#======================#
sub usage {
### Returns usage information for current script in a string
<<USAGE;
Usage: perl $0 [OPTIONS] <FASTA>

Description: Computes nucleotide counts, sequence length and GC content from one or more FASTA files.

Options:
	--trim-id               From each ID, remove everything following the first whitespace (default: do not trim).
	--usage | --help        Show this screen and exit.
	--version               Show version information and exit.
	--verbose               Print log information to STDERR.

Comments:
	- Only marginal validation of FASTA format performed. In particular, empty sequences are not handled.
	- U and T are both counted. GC content is calculated as the fraction of (G + C) / (A + C + G + T + U) * 100.
	- Specified length values equal the number of letters in the sequence minus whitespace.

Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Sep 22, 2014.
Version 1.0.1 (May 20, 2015)
USAGE
}
#-----------------------#
sub version {
### Returns version information for current script in a string
<<VERSION;
$0 version 1.0.1 (May 20, 2015)
Written by Alexander Kanitz, Biozentrum, Univerity of Basel on Sep 22, 2014.
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
my $verbose = 0;

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
my @seq = ();
my $seq = "";
my @stats = ();

	#---> BODY <---#

	#---> Print header line <---#
	print STDOUT "ID\tPercent_GC\tLength\tCount_ACGTU\tCount_A\tCount_C\tCount_G\tCount_T\tCount_U\n";
	
	#---> Loop over input lines <---#
	while (my $line = <>) {

		#---> Remove trailing newline character <---#
		chomp $line;
	
		#---> IF line is header line... <---#	
		if ($line =~ /^>/) {
	
			#---> Test IF there is already a sequence (should be TRUE for all but the first entry) <---#
			if ( scalar @seq ) {
	
					#---> Join all entries in the sequences array and empty array <---#
					$seq = join "", @seq;
					@seq = ();
	
					#---> Compute nucleotide stats from joined sequence <---#
					@stats = @{&nucleotide_stats($seq)};
					
					#---> Print nucleotide stats <---#
					print STDOUT (join "\t", @stats), "\n";
				
			}
	
			#---> Trim identifier (remove trailing ">" and - optionally - everything after the first whitespace) <---#
			if ($trim_id) {
				$line =~ s/^>(.+?)(\s.+)?$/$1/g;
			}
			else {
				$line =~ s/^>//;
			}
	
			#---> Print out ID
			print STDOUT "$line\t";

		}
	
		#---> ...ELSE push line to sequence array <---#
		else { push @seq, $line };

	}

	#---> Complete last entry: Join all entries in the sequences array <---#
	$seq = join "", @seq;
	
	#---> Complete last entry: Compute nucleotide stats from joined sequence <---#
	@stats = @{&nucleotide_stats($seq)};
	
	#---> Complete last entry: Print nucleotide stats <---#
	print STDOUT (join "\t", @stats), "\n";

#---> STATUS MESSAGE <---#
print STDERR "Done.\n" if $verbose;

#---> PROGRAM EXIT <---#
exit 0;
#===========#
#  // MAIN  #
#===========#


#==================#
#  SUBROUTINES //  #
#==================#
sub nucleotide_stats {
## Description: Counts nucleotides and computes GC content from a nucleotide sequence.
## Accepts: A string containing a nucleotide sequence.
## Returns: A reference to an array of the 1. GC content, 2. sequence length, 3. total number of valid nucleobases, 4.-8. counts of nucleobases A, C, G, T & U
## Dependencies: n/a
## Type: Generic

	#---> PASS ARGUMENTS ---#
	my $seq = shift();
	my @return = ();

	#---> SUBROUTINE VARIABLES <---#
	my $gc_content = "NA";

	#---> BODY <---#

		#---> Remove whitespace from sequence <---#
		$seq =~ s/\s//g;

		#---> Count number of occurrences of each valid nucleotide base in sequence <---#
		my $a_count = () = $seq =~ m/a/ig;
		my $c_count = () = $seq =~ m/c/ig;
		my $g_count = () = $seq =~ m/g/ig;
		my $t_count = () = $seq =~ m/t/ig;
		my $u_count = () = $seq =~ m/u/ig;

		#---> Get sequence length <---#
		my $length = length($seq);

		#---> Get total number of valid nucleotides <---#
		my $acgtu_count = $a_count + $c_count+ $g_count + $t_count + $u_count;

		#---> Compute GC content <---#
		$gc_content = (( $c_count + $g_count) / $acgtu_count * 100 ) if $acgtu_count > 0;

		#---> Assemble return string <---#		
		push @return, $gc_content, $length, $acgtu_count, $a_count, $c_count, $g_count, $t_count, $u_count;

	#---> RETURN VALUE <---#
	return \@return;

}
#==================#
#  // SUBROUTINES  #
#==================#
