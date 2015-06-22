#!/import/bc2/soft/bin/perl5/perl -w
#!/import/bc2/soft/bin/perl5/perl -w
#===============================================================================
#       USAGE: removeAdaptor.pl seq-copies-file alignments
#       PURPOSE: given alignment of reads to adaptor, the script removes adaptor fragments from the reads
#       REQUIREMENTS:  Perl 5.8.1
#       COMMENTS:
#       AUTHOR:  Philipp Berninger, Chris Rodak
#       COMPANY:  Biozentrum, University of Basel
#===============================================================================

use strict;
use warnings;

if ( @ARGV != 1 ) {
  die "usage: ./removeAdaptor.pl alignments\n";
}

my %newseqs;

#seq     TGAGGTAGTAGGTTGTATGGTTTCGTA----------------
#ada     ----------------------TCGTATGCCGTCTTCTGCTTG
open( IN, $ARGV[0] ) || die "cannot open $ARGV[0]!\n";
while (<IN>) {
  chomp($_);
  if ( $_ =~ /^seq\S*\s+(\S+)/ ) {
    my @d = split( /\s+/, $_ );
    my $s1 = $1;
    $_ = <IN>;
    if ( $_ =~ /^\S+\s+(\S+)/ ) {
      my $s2      = $1;
      my $adaptor = $s2;
      $adaptor =~ s/[\s\-]//g;
      my $adaptorlength = length($adaptor);
      my @c1            = split( //, $s1 );
      my @c2            = split( //, $s2 );
      my $refseq        = $s1;
      $refseq =~ s/\-//g;
      my $i = 0;
      my ( $match, $mismatch ) = ( 0, 0 );

      while ( $i < @c2 && $c2[$i] =~ /\-/ ) {
        $i++;
      }
      my $j = $#c1;
      while ( $j >= 0 && $c1[$j] =~ /\-/ ) {
        $j--;
      }
      for ( my $k = $i ; $k <= $j ; $k++ ) {
        if ( $c1[$k] eq $c2[$k] ) {
          $match++;
        } else {
          $mismatch++;
        }
      }
      if ( ( $match - $mismatch ) >= 3 && $mismatch <= 1 ) {

        #this is basically the recognizable adaptor case
        my $extractedseq = substr( $s1, 0, $i );
        $extractedseq =~ s/\-//g;

        #chop off all the Ns from the end of the Sequence
        if ( $extractedseq =~ /(N*)$/ ) {
          $extractedseq = substr( $extractedseq, 0, length($extractedseq) - length($1) );
        }

        #chop off all the Ns from the begining of the Sequence
        if ( $extractedseq =~ /^(N*)/ ) {
          $extractedseq = substr( $extractedseq, length($1), length($extractedseq) - length($1) );
        }
        if ($extractedseq) {
          print $extractedseq . "\t" . $d[2] . "\n";
        }
      } elsif ( ( $match == 2 || $match == 1 ) && $mismatch < 1 ) {

        #no mismatch, but only 1-2 nucleotides matching adaptor
        my $extractedseq = substr( $s1, 0, $i );
        $extractedseq =~ s/\-//g;

        #chop off all the Ns from the end of the Sequence
        if ( $extractedseq =~ /(N*)$/ ) {
          $extractedseq = substr( $extractedseq, 0, length($extractedseq) - length($1) );
        }

        #chop off all the Ns from the begining of the Sequence
        if ( $extractedseq =~ /^(N*)/ ) {
          $extractedseq = substr( $extractedseq, length($1), length($extractedseq) - length($1) );
        }
        if ($extractedseq) {
          print $extractedseq . "\t" . $d[2] . "\n";
        }
      } else {

        #check if there is extra sequence after the adaptor
        #deals with some cases that we observed at some point
        $j = $#c2;
        while ( $j >= 0 && $c2[$j] =~ /\-/ ) {
          $j--;
        }
        for ( my $k = $i ; $k <= $j ; $k++ ) {
          if ( $c1[$k] eq $c2[$k] ) {
            $match++;
          } else {
            $mismatch++;
          }
        }
        if ( $match / $adaptorlength >= 0.85 ) {
          my $extractedseq = substr( $s1, 0, $i );
          $extractedseq =~ s/\-//g;

          #chop off all the Ns from the end of the Sequence
          if ( $extractedseq =~ /(N*)$/ ) {
            $extractedseq = substr( $extractedseq, 0, length($extractedseq) - length($1) );
          }

          #chop off all the Ns from the begining of the Sequence
          if ( $extractedseq =~ /^(N*)/ ) {
            $extractedseq = substr( $extractedseq, length($1), length($extractedseq) - length($1) );
          }
          if ($extractedseq) {
            print $extractedseq . "\t" . $d[2] . "\n";
          }
        } else {

          #read without recognizable adaptor or adaptor in some other configuration than
          #RRRRRRRRRRRDDDDDDD or RRRRRRRRRRRRDDDDDDDDJJJJJJJ (R - RNA, D - adaptor, J - junk)
          #chop off all the Ns from the end of the Sequence
          if ( $refseq =~ /(N*)$/ ) {
            $refseq = substr( $refseq, 0, length($refseq) - length($1) );
          }

          #chop off all the Ns from the begining of the Sequence
          if ( $refseq =~ /^(N*)/ ) {
            $refseq = substr( $refseq, length($1), length($refseq) - length($1) );
          }
          print $refseq . "\t" . $d[2] . "\n";
        }
      }
    }
  }
}
close IN;
