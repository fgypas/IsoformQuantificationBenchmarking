#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: sam_add_mm_linked_list_encoding.pl
### Created: Jan 26, 2015
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Requirements: GetOpt::Long
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
my $usage = 0;
my $verbose = 0;
my $in = '';
my $out = '';
my $no_head = 0;
my $link_mm = 0;
my $nh = 0;
my $hi = 0;
my $mapq = 0;
my $qual = 0;
my $qual_char = "I";
my $remove_custom = 0;
my $options_result = GetOptions (
        'usage|help' => \$usage,
        'verbose' => \$verbose,
        'discard-header' => \$no_head,
	'link-mm' => \$link_mm,
	'nh' => \$nh,
	'mapq' => \$mapq,
	'qual' => \$qual,
	'qual-char=s' => \$qual_char,
	'remove-custom-tags' => \$remove_custom,
	'enforce-hi' => \$hi,
        'in=s' => \$in,
        'out=s' => \$out
);

## Die if command line parsing was not successful or required arguments are missing
die $usage_info if $usage || !$options_result;

## Die if specified argument to --qual is illegal
die "[ERROR] Specified argument to --qual is illegal. Only a single ASCII character is allowed.\n$usage_info" if length $qual_char > 1;

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print STDERR "Starting '$0'...\n" if $verbose;

#---> BODY <---#

	## Call main subroutine
	&add_tags_sam($in, $out, $link_mm, $nh, $hi, $mapq, $qual, $qual_char, $remove_custom, $no_head);

#---> STATUS MESSAGE <---#
print STDERR "Done.\n" if $verbose;

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
'Usage: perl ./sam_filter_gen.pl [OPTIONS] --in [SAM] --out [SAM]

Description: Adds tags allowing linked-list encoding of reads that align to more than one reference
             sequence ("RNAME") to each entry of a SAM file sorted by QNAME. The following tags can
             be added: CC/CP (RNAME and POS of next alignment of the same read/QNAME), NH (total
             number of hits/alignments of a given read/QNAME), HI (hit index of the current alignment
             out of all alignments for a given read/QNAME; 0-based). Additionally, the MAPQ-field can
             be updated in a TopHat-like fashion and the QUAL field replaced with a "dummy" Phred
             score. Existing fields with the same identifiers will be replaced. Refer to the SAM
             format specifications (currently at http://samtools.github.io/hts-specs/SAMv1.pdf) for
             further details.

Arguments:
  --in [SAM]            Input SAM file sorted by QNAME (default: STDIN).
  --out [SAM]           Output SAM file (default: STDOUT).
  --discard-header      Discard SAM header (only print alignments).
  --link-mm             If more than one alignment is reported for a given read (QNAME), add the
                        RNAME and POS field values R and P of the next alignment with the same QNAME
                        field value as TAG fields of the form "CC:Z:R" and "CP:i:P", respectively.
                        Note that if the next alignment has the same RNAME as the current one, R will
                        be set to "=". Moreover, a 0-based index I of the current alignment out of
                        all alignments with the same QNAME field value as a TAG field of the form
                        "HI:i:I". The first out of a set of multiple alignments sharing the same
                        QNAME is chosen as the beginning of the resulting "linked list" and will
                        receive index I = 0. No CC or CP TAG fields will be added to the last
                        alignment of each such set. For all but the first alignment, the "secondary
                        alignment" FLAG bit 0x100 is set (256 is added to current decimal FLAG field
                        value).
  --nh                  Add the total number N of all alignments that share the same QNAME field
                        value as a TAG field of the form "NH:i:N".
  --enforce-hi          Add the HI TAG field even if --link-mm is not set or there are no secondary
                        alignments.
  --mapq                Add a Phred-score-like value to the MAPQ field that is derived from the
                        number of alignments reported for a given read in the manner as is used by
                        the TopHat program. One of four values is used for all alignments of a given
                        read:
                        - 50 (only one alignment)
                        - 3 (two alignments)
                        - 1 (three or four alignments)
                        - 0 (five or more alignments)
  --qual                Substitute the QUAL field value with a "dummy" Phred-score-based sequencing
                        quality string. It consists of n repeats of the ASCII character specified via
                        --qual-char, where n is the length of the sequence in the SEQ field. 
  --qual-char [CHAR]	When --qual is specified, use the specified ASCII character to generate the 
                        "dummy" sequencing quality string (default: "I", corresponding to the highest
                        quality score in the Sanger FASTQ format). This option is ignored if --qual
                        is not specified.
  --remove-custom-tags  Remove all custom TAG fields. According to the SAM specification, these
                        include all whose two-character identifier string starts with X, Y, Z or
                        contains a lower case character at either position.
  --usage|help          Show this information and die.
  --verbose             Print status messages to STDERR.

Notes:
CAUTION: The correct sorting order is not validated! To ensure the correct sorting is maintained,
         execute "samtools sort <IN> <OUT_PREFIX>" before running this script.
CAUTION: Only marginal validation of the input file type/format order performed!

Version 1.0 (2015-01-26)
Written by Alexander Kanitz on 2015-01-26
';
}
#-----------------------#
sub add_tags_sam {
### Function: Adds tags allowing linked-list encoding of reads that align to more than one locus ("multimappers") to each entry of a SAM file sorted by QNAME. The following tags can be added: CC/CP (RNAME and POS of next alignment of the same read/QNAME), NH (total number of hits/alignments of a given read/QNAME), HI (hit index of the current alignment out of all alignments for a given read/QNAME; 0-based). Existing TAG fields with the same identifiers will be replaced. Refer to the SAM format specifications (currently at http://samtools.github.io/hts-specs/SAMv1.pdf) for further details.
### Accepts: 1. Input file [FILE|SAM]; 2. Output file [FILE|SAM]; 3. CC/CP switch: 0 -> inactive, 1 -> CC and CP TAG fields will be added; 4. NH switch: 0 -> inactive, 1 -> NH TAG field will be added; 5. HI switch: 0 -> inactive, 1 -> HI TAG field will be added if more than one alignment in array of hashes, 2 -> HI TAG will always be added; 6. Header switch: FALSE = do not print header, TRUE = print header; 7. Header file (prepends SAM records in output)
### Returns: n/a
### Dependencies: n/a
### Type: Specialized
        #---> PASS ARGUMENTS ---#
        my ($in, $out, $link_mm, $nh, $hi, $mapq, $qual, $qual_char, $remove_custom, $no_head) = @_;

        #---> STATUS MESSAGE <---#
        print STDERR "Adding/updating TAG fields to SAM file '$in'..." . "\n" if $verbose;

        #---> SUBROUTINE VARIABLES <---#
        my $regex_header = '^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$';
        my $regex_comment = '^\@CO.*$';
        my $last_line;                                                                                                                                                                                                          # holds the last line
        my $last_id;                                                                                                                                                                                                            # holds the last distinct QNAME/read ID
        my @AoH;                                                                                                                                                                                                                        # holds references to hashes containing the field values of the last lines that share the same QNAME/read ID
        my @field_keys = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL/;                                         # holds the bareword keys for the hashes containing field values
        my $final_record = 0;

        #---> BODY <---#

                #---> Open input and output filehandles <---#
		my ($in_fh, $out_fh);
		my $in_open = open $in_fh , '<', $in or die "[ERROR] Could not open file '$in'!\n" if $in;
		open $in_fh, '<&STDIN' unless $in_open;
		my $out_open = open $out_fh, '>', $out or die "[ERROR] Could not open file '$out'!\n" if $out;
		$out_fh = *STDOUT unless $out_open;

                #---> Traverse non-header lines <---#
                while (my $line = <$in_fh>) {

			#---> Print header <---#
                        if ( $line =~ /$regex_header/ || $line =~ /$regex_comment/ ) {
                                print $out_fh $line unless $no_head;
				next;
                        }
                        
			#---> Field values to hash <---#
			chomp $line;
                        my @field_values = split "\t", $line;
                        die "[ERROR] Input file does not look like a valid SAM file!" unless scalar @field_values >= 11;        # Assert presence of at least 11 fields
                        my %fields;
                        @fields{@field_keys} = @field_values[0 .. 10];
                        foreach (@field_values[11 .. $#field_values]) {
                                my ($tag, $value) = split ":", $_, 2;
                                $fields{$tag} = $value;
                        }

	                #---> Update QUAL field <---#
                	$fields{"QUAL"} = $qual_char x length $fields{"SEQ"} if $qual;

	                #---> Remove custom tags <---#
			if ($remove_custom) {
				my @custom_keys = grep { /^[XYZa-z][A-Za-z0-9]$/ } keys %fields;
				delete $fields{$_} for @custom_keys;
			}

                        #---> Manage AoH: Grow if QNAMEs identical, else compare AoH entries, print record(s) and reset AoH <---#
                        if ( defined $last_id && $fields{"QNAME"} ne $last_id ) {
                                die "[ERROR] SAM file appears to be corrupt!" unless scalar @AoH > 0;                           # Assert integrity of SAM records
                                @AoH = @{&sam_add_mm_encodings(\@AoH, $link_mm, $nh, $hi, $mapq)};
                                my @out_lines = @{&sam_AoH_join_records(\@AoH)};
				print $out_fh $_ foreach @out_lines;                                                            # Print entries
                                @AoH = ();                                                                                      # Clear array
                        }
                        push @AoH, \%fields;
                        $last_id = $fields{"QNAME"};
                        $last_line = $line;

                }

                #---> Account for final record(s) separately due to EOF <---#
                die "[ERROR] SAM file appears to be corrupt or empty!" unless scalar @AoH > 0;                                  # Assert integrity of SAM records
                @AoH = @{&sam_add_mm_encodings(\@AoH, $link_mm, $nh, $hi, $mapq)} if scalar @AoH > 1;
                my @out_lines = @{&sam_AoH_join_records(\@AoH)};
		print $out_fh $_ foreach @out_lines;                                                                            # Print entries

                #---> Close input and output filehandles <---#
                close $out_fh if $out_open;
                close $in_fh if $in_open;

        #---> END BODY <---#

        #---> STATUS MESSAGE <---#
        print STDERR "Written modified SAM file to '$out'.\n" if $verbose;

        #---> RETURN VALUE <---#
        return 0;
}
#-----------------------#
sub sam_add_mm_encodings {
### Function: Adds tags allowing linked-list encoding of reads that align to more than one locus ("multimappers") to each entry of a SAM file sorted by QNAME. The following tags can be added: CC/CP (RNAME and POS of next alignment of the same read/QNAME), NH (total number of hits/alignments of a given read/QNAME), HI (hit index of the current alignment out of all alignments for a given read/QNAME; 0-based). Existing TAG fields with the same identifiers will be replaced. Refer to the SAM format specifications (currently at http://samtools.github.io/hts-specs/SAMv1.pdf) for further details.
### Accepts: 1. Array of hashes of SAM records with identical QNAME field (as generated by the subroutine 'add_tags_sam', written by Alexander Kanitz, 26-JAN-2015); 2. CC/CP switch: 0 -> inactive, 1 -> CC and CP TAG fields will be added; 3. NH switch: 0 -> inactive, 1 -> NH TAG field will be added; 4. HI switch: 0 -> inactive, 1 -> HI TAG field will be added if more than one alignment in array of hashes, 2 -> HI TAG will always be added
### Returns: Reference to array of hashes
### Dependencies: n/a
### Type: Generic
        #---> PASS ARGUMENTS ---#
        my ($AoH_ref, $link_mm, $nh, $hi, $mapq) = @_;

        #---> BODY <---#

                #---> Add CC, CP, HI TAG & SET SECONDARY ALIGNMENT BIT <---#
		if ( $link_mm && scalar @$AoH_ref > 1 ) {
			for my $i (0 .. $#{$AoH_ref} ) {					
				$$AoH_ref[$i]->{"FLAG"} = $$AoH_ref[$i]->{"FLAG"} + 256 unless $i == 0;
				$$AoH_ref[$i]->{"HI"} = "i:" . $i;
				$$AoH_ref[$i]->{"CC"} = $$AoH_ref[$i+1]->{"RNAME"} eq $$AoH_ref[$i]->{"RNAME"} ? "Z:=" : "Z:" . $$AoH_ref[$i+1]->{"RNAME"} unless $i == $#{$AoH_ref};
				$$AoH_ref[$i]->{"CP"} = "i:" . $$AoH_ref[$i+1]->{"POS"} unless $i == $#{$AoH_ref};
			}
		}

                #---> Add HI TAG <---#
		if ($hi) {
			for my $i (0 .. $#{$AoH_ref} ) {					
				$$AoH_ref[$i]->{"HI"} = "i:" . $i;
			}
		}

                #---> Add NH TAG <---#
		if ($nh) {
			foreach my $hash_ref (@$AoH_ref) {
				$hash_ref->{"NH"} = "i:" . scalar @$AoH_ref;			
			}
		}

		#---> Update MAPQ field <---#
		if ($mapq) {
                        foreach my $hash_ref (@$AoH_ref) {
                                $hash_ref->{"MAPQ"} = 50 if scalar @$AoH_ref == 1;
				$hash_ref->{"MAPQ"} = 3 if scalar @$AoH_ref == 2;
				$hash_ref->{"MAPQ"} = 1 if scalar @$AoH_ref == 3 || scalar @$AoH_ref == 4;
				$hash_ref->{"MAPQ"} = 0 if scalar @$AoH_ref >= 5;
                        }
		}

        #---> RETURN VALUE <---#
        return $AoH_ref;
}
#-----------------------#
sub sam_AoH_join_records {
### Function: Joins the fields of SAM records stored in an array of hashes in the proper order (additional tags in alphanumerical order)
### Accepts: Array of hashes of SAM records with identical QNAME field (as generated by the subroutine 'filter_sam', written by Alexander Kanitz, 29-AUG-2013)
### Returns: Array of strings, one element for each record, appended with a newline character for easy printing
### Dependencies: n/a
### Type: Generic
        #---> PASS ARGUMENTS ---#
        my $AoH_ref = shift;

        #---> SUBROUTINE VARIABLES <---#
        my @field_keys = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL/;                        # holds the bareword keys for the hashes containing field values in the right order
        my @records;

        #---> BODY <---#

                foreach my $record (@$AoH_ref) {
                        my @fields_ordered;
                        foreach my $key (@field_keys) {
                                push @fields_ordered, $record->{$key};
                                delete $record->{$key};
                        }
                        foreach my $extra_field (sort keys %$record) {
                                push @fields_ordered, ( $extra_field . ":" . $record->{$extra_field} );
                        }
                        push @records, ( join( "\t", @fields_ordered) . "\n" );
                }

        #---> RETURN VALUE <---#
        return \@records;
}
#=======================#
#    SUBROUTINES END    #
#=======================#

