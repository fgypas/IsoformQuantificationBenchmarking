use strict;
use warnings;
use Getopt::Long;

my $help;
my $valid;
my $debug;
my $strict;
my $exclude;
my $correction = 0;
GetOptions(
  "strict"       => \$strict,
  "exclude=s"    => \$exclude,
  "correction=i" => \$correction,
  "valid=s"      => \$valid,
  "debug"        => \$debug,
  "help"         => \$help,
  "h"            => \$help
);

my $showhelp = 0;
if ( defined $strict ) {
  $showhelp = 1 if ( not defined $valid );
}
$showhelp = 1 if ( defined $help );

if ($showhelp) {
  print "Usage: zcat file | perl $0 --valid=file.ids.gz\n\n";
  exit 1;
}

my %V = ();
if ($strict) {
  if ( not -e $valid ) {
    print STDERR "[ERROR] File '$valid' not found.\n";
    exit;
  }

  if ( $valid !~ m/.gz/ ) {
    print STDERR "[ERROR] File '$valid' does not look like a gzipped file.\n";
    exit;
  }
  open( IN, "gunzip -c $valid |" );
  while (<IN>) {
    chomp;
    $V{$_} = 1;
  }
  close(IN);
}

# %ordered stores an ordered set of chromosome names
# starting with chr1 ... chrM
my %ordered     = ();
my $idx_ordered = 0;
foreach my $strand ( '+', '-' ) {
  foreach my $i ( 1 .. 100, 'X', 'Y' ) {
    $idx_ordered++;
    $ordered{"chr$i:$strand"} = $idx_ordered;
    $idx_ordered++;
    $ordered{"chr$i\_random:$strand"} = $idx_ordered;
  }
}

my %hash      = ();
my %hashdebug = ();
while (<STDIN>) {
  chomp;
  my @F = split(/\t/);
  ( my $readname, my $readcount ) = split( /-/, $F[3] );
  next if ( $F[0] eq $exclude );
  if ($strict) {
    next if ( not defined $V{ $F[3] } );
    my $alignflag = 0;
    if ( $F[5] eq '+' ) {
      if ( $F[6] =~ m/^\d+$/ ) {

        # perfect match: okay
        $alignflag = 1;
      } elsif ( $F[6] !~ m/\d+$/ ) {

        # end badly aligned: skip
      } else {
        if ( $F[6] =~ m/(.*\D)(\d+)$/ ) {
          if ( $2 > 3 ) {
            $alignflag = 1;

            #print STDERR "$_\n";
          }
        }
      }
    } else {
      if ( $F[6] =~ m/^\d+$/ ) {

        # perfect match: okay
        $alignflag = 1;
      } elsif ( $F[6] !~ m/^\d+/ ) {

        # end badly aligned: skip
      } else {
        if ( $F[6] =~ m/^(\d+)(\D.*)/ ) {
          if ( $1 > 3 ) {
            $alignflag = 1;

            #print STDERR "$_\n";
          }
        }
      }
    }
    next if ( $alignflag == 0 );
  }

  my $key = "$F[0]:$F[5]";
  $F[1]++;

  # add the chromosome/strand key if not
  # already in %ordered
  if ( not defined $ordered{$key} ) {
    $idx_ordered++;
    $ordered{$key} = $idx_ordered;
  }

  $hash{$key} = {} if ( not defined $hash{$key} );
  if ($debug) {
    $hashdebug{$key} = {} if ( not defined $hashdebug{$key} );
  }
  if ( $F[5] eq '+' ) {
    $hash{$key}->{ $F[2] } += $readcount;
    if ($debug) {
      $hashdebug{$key}->{ $F[2] } = [] if ( not defined $hashdebug{$key}->{ $F[2] } );
      push @{ $hashdebug{$key}->{ $F[2] } }, $F[3];
    }
  } else {
    $hash{$key}->{ $F[1] } += $readcount;
    if ($debug) {
      $hashdebug{$key}->{ $F[1] } = [] if ( not defined $hashdebug{$key}->{ $F[1] } );
      push @{ $hashdebug{$key}->{ $F[1] } }, $F[3];
    }
  }
}

my $idx = 0;
foreach my $key ( sort { $ordered{$a} <=> $ordered{$b} } keys %ordered ) {
  next if ( not defined $hash{$key} );
  my ( $chr, $strand ) = split( /:/, $key );
  foreach my $site ( sort { $a <=> $b } keys %{ $hash{$key} } ) {
    my $start = $site + $correction - 1;
    my $end   = $site + $correction;
    if ( $strand eq '-' ) {
      $start = $site - $correction - 1;
      $end   = $site - $correction;
    }
    my $score = $hash{$key}->{$site};
    $idx++;
    my $additional = '';
    if ($debug) {
      $additional = "\t" . join( ",", @{ $hashdebug{$key}->{$site} } );
    }
    print "$chr\t$start\t$end\t$idx\t$score\t$strand$additional\n";
  }
}

