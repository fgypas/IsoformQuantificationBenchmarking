use strict;
use warnings;

my %READCOUNT = ();
open( F, "gunzip -c $ARGV[0] |" );
while (<F>) {
  chomp;
  my @F = split(/\t/);
  $READCOUNT{ $F[3] }++;
}
close(F);

open( F, "gunzip -c $ARGV[0] |" );
while (<F>) {
  chomp;
  my @F = split(/\t/);
  ( my $readname, my $readcount ) = split( /-/, $F[3] );
  $readcount = sprintf( "%.2f", $readcount / $READCOUNT{ $F[3] } );
  my $pos = ( $F[5] eq '+' ) ? $F[2] : $F[1] + 1;
  my $start = $pos - 1;
  print "$F[0]\t$start\t$pos\tM$readname\t$readcount\t$F[5]\n";
}
close(F);

