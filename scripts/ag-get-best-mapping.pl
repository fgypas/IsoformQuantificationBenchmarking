use strict;
use warnings;

my %hash = ();

foreach my $i ( 0 .. $#ARGV ) {
  if ( not -e $ARGV[$i] ) {
    print STDERR "[ERROR] File '$ARGV[$i]' not found.\n";
    exit 1;
  }
  print STDERR "[INFO] Reading file '$ARGV[$i]'.\n";
  open( F, "gunzip -c $ARGV[$i] |" );
  while (<F>) {
    chomp;
    my @F = split( /\t/, $_ );
    my $pos = "$F[0]:$F[1]:$F[2]:$F[5]";
    if ( defined $hash{ $F[3] } ) {

      # lets check if the exact same location is already there
      my $flag_present      = 0;
      my $min_edit_distance = 10e6;
      foreach my $entry ( @{ $hash{ $F[3] } } ) {
        my $postmp = "$entry->[0]:$entry->[1]:$entry->[2]:$entry->[5]";
        $flag_present = 1 if ( $postmp eq $pos );
        $min_edit_distance = $entry->[4] if ( $entry->[4] < $min_edit_distance );
      }
      next if ( $flag_present == 1 );
      if ( $F[4] == $min_edit_distance ) {
        push @{ $hash{ $F[3] } }, \@F;
        next;
      }

      if ( $F[4] < $min_edit_distance ) {

        #print STDERR "$i: $F[4] < $min_edit_distance\n";
        $hash{ $F[3] } = [];
        push @{ $hash{ $F[3] } }, \@F;
      }
    } else {
      $hash{ $F[3] } = [];
      push @{ $hash{ $F[3] } }, \@F;
    }
  }
  close(F);
}

foreach my $key ( keys %hash ) {
  # only output reads that are unique, either on the genome or on the transcriptome
  # discard reads that map to genome and transcriptome with the same edit-distance
  next if ( scalar @{ $hash{$key} } > 1);
  foreach my $entry ( @{ $hash{$key} } ) {
    print join( "\t", @{$entry} ), "\n";
  }
}
