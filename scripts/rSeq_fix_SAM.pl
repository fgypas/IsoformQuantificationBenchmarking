################################################################################
## Author: Foivos Gypas							      ##
## Date: 07-OCT-2014							      ##
## Zavolan Group, Biozentrum, University of Basel, Switzerland		      ##
################################################################################

use strict;
use warnings;

my $dictionary_path = $ARGV[0]; 
my $sam_path = $ARGV[1];

# Create hash TranscriptID ==> GeneID$$TranscriptID
my %dict;
open(IN, $dictionary_path) or die "Can't open $!";
while(<IN>){
	my $line = $_; # read line by line
	chomp $line; # remove \n
	my @f=split('\t', $line);
	$dict{$f[0]}=$f[1]; #fill dictionary
}
close IN;

# Read SAM file and rename TranscriptID to GeneID$$TranscriptID
open(IN, $sam_path) or die "Can't open $!";
while(<IN>){
	# Rename header in SAM file
	if($_ =~ m/^\@SQ\t/){
		my @sp_line = split("\t",$_);
		my $l = substr $sp_line[1], 3, 15;
		$sp_line[1] =  $dict{$l};	
		$sp_line[1] = "SN:" . $sp_line[1];
		print join("\t", @sp_line);
	}
	elsif($_ =~ m/^\@/){
		print $_;
	}
	# Rename body
	else{
		my @sp_line = split("\t",$_);
		my $l = $sp_line[2];
		$sp_line[2] =  $dict{$l};
		print join("\t", @sp_line);
	}
}
close IN;
