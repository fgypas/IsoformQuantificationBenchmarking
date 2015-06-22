use strict;
use warnings;

my %hash = ();
open( F, $ARGV[0] );
while (<F>) {
  chomp;
  if ( $_ =~ m/^clusterid/ ) {
    print "$_\n";
    next;
  }
  my @F = split(/\t/);
  $hash{ $F[0] } = [ 0, 0 ] if ( not defined $hash{ $F[0] } );
  $hash{ $F[0] }->[0] += $F[1];
  $hash{ $F[0] }->[1] += $F[2];
}
close(F);

open( F, $ARGV[1] );
while (<F>) {
  chomp;
  if ( $_ =~ m/^clusterid/ ) {
    next;
  }
  my @F = split(/\t/);
  $hash{ $F[0] } = [ 0, 0 ] if ( not defined $hash{ $F[0] } );
  $hash{ $F[0] }->[0] += $F[1];
  $hash{ $F[0] }->[1] += $F[2];
}
close(F);

foreach my $key ( sort keys %hash ) {
  my @X = @{ $hash{$key} };
  foreach my $idx ( 0 .. $#X ) {
    $X[$idx] = sprintf( "%.2f", $X[$idx] );
  }
  print "$key\t", join( "\t", @X ), "\n";
}
