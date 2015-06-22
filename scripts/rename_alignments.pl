################################################################################
## Author: Foivos Gypas							      ##
## Date: 08-OCT-2014							      ##
## Zavolan Group, Biozentrum, University of Basel, Switzerland		      ##
################################################################################

use strict;
use warnings;

#my $dictionary_path = $ARGV[0];
my $sam_path = $ARGV[0];

# Initialize
my $current_key = '';
my $current_value ='';
my @chars = ("A".."Z", "a".."z", "0".."9" );
my @tmp = ();

# Read SAM file and rename
open(IN, $sam_path) or die "Can't open $!";
while(my $line = <IN>){

	# Rename header in SAM file.
	if($line =~ m/^\@/){
		print $line;
	}
	# Rename body
	else{
		
		my @sp_line = split("\t",$line); # Split line
		my $name = $sp_line[0]; # Get old read name

		# If same key exists
		if($name eq $current_key){
 			$sp_line[0] = $current_value; #$tmp_dict{$name};
	                print join("\t", @sp_line);
		}
		else{	
			# Generate new read name
			#my $r = `~/bin/mktemp -u XXXXXXXXXXXXXXXXXXXX`;
			#chomp $r;
			#my $new_name = join("","s", $r); # add an s a the beginning of the read 
			@tmp = ();
			push(@tmp, 's');
			push(@tmp, $chars[rand @chars]) for 1..20;
			my $new_name = join('', @tmp);


			$current_key = $name;
			$current_value = $new_name;

			$sp_line[0] = $current_value;
			print join("\t", @sp_line);
		}
	}
}
close IN;
