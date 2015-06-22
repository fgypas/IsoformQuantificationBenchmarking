use strict;
use warnings;
use Getopt::Long;

my @files = ();
my @names = ();
my $help;
my $tpm;
my $pipe;
my $clusterfile;
my $minLevel = 2;
my $scoreAsReadCount;
my $savedRejected;
my $chrM;
my $IPfile;

GetOptions(
  "sample=s"         => \@files,
  "clusters=s"       => \$clusterfile,
  "scoreAsReadCount" => \$scoreAsReadCount,
  "savedRejected"    => \$savedRejected,
  "minLevel=i"       => \$minLevel,
  "tpm"              => \$tpm,
  "IPfile:s"         => \$IPfile,
  "chrM:s"           => \$chrM,
  "name=s"           => \@names,
  "pipe"             => \$pipe,
  "help"             => \$help,
  "h"                => \$help
);

# process options
my $showhelp = 0;

if ( !$pipe ) {
  $showhelp = 1 if ( $#files == -1 );
  $showhelp = 1 if ($help);
}

if ( $pipe and $savedRejected ) {
  print STDERR "[ERROR] --savedRejected not possible with --pipe.\n";
  exit;
}

$showhelp = 1 if ( not defined $clusterfile );

if ($showhelp) {
  print STDERR "Usage: $0 --clusters=file.bed --sample=file1.bed --name=name1 > out.tsv\n\n";
  print STDERR "Mandatory parameters\n\n";
  print STDERR "   --clusters=file      BED file obtained from polyasite.unibas.ch\n\n";
  print STDERR "   --sample=file        BED file (multiple values allowed)\n\n";
  print STDERR "\nOptional parameters\n\n";
  print STDERR "   --name=text          name of the sample (multiple values allowed)\n\n";
  print STDERR "   --scoreAsReadCount   Use the 5th field in the sample BED file as read count\n\n";
  print STDERR "   --minLevel=2         Minimal support level for clusters\n\n";
  print STDERR "   --savedRejected      Write rejected reads to a BED file\n\n";
  print STDERR "   --IPfile=file        Use the PolyASite IP annotation BED file\n\n";
  print STDERR
    "   --chrM=chrM          Do not count read originating from the mitochondrial genome\n\n";
  print STDERR "   --tpm                Normalize reads to tags per million (TPM)\n\n";
  print STDERR "   --pipe               Read from STDIN instead of a file\n\n";
  exit 1;
}

if ( $#names > -1 and $#files != $#names ) {
  print STDERR "[ERROR] # names does not match # samples.\n";
  exit 1;
}

if ( not -e $clusterfile ) {
  print STDERR "[ERROR] File '$clusterfile' does not exist.\n";
  exit 1;
}
if ( -d $clusterfile ) {
  print STDERR "[ERROR] '$clusterfile' is not a file.\n";
  exit 1;
}

foreach my $file (@files) {
  if ( not -e $file ) {
    print STDERR "[ERROR] File '$file' does not exist.\n";
    exit 1;
  }
  if ( -d $file ) {
    print STDERR "[ERROR] '$file' is not a file.\n";
    exit 1;
  }
}

# check if gunzip is installed
my $gunzip    = `which gunzip 2>&1`;
my $no_gunzip = 0;
if ( $gunzip =~ m/no\sgunzip\sin/ ) {
  $no_gunzip = 1;
  print STDERR "[INFO] Unzipping of .gz files not possible.\n";
}

foreach my $file (@files) {
  if ( $file =~ m/\.gz$/ and $no_gunzip == 1 ) {
    print STDERR "[INFO] Unzipping is not possible.\n";
    exit 1;
  }
}

my %IP = ();
if ($IPfile) {
  if ( not -e $IPfile ) {
    print STDERR "[ERROR] File '$IPfile' does not exist.\n";
    exit 1;
  }
  my $fh;
  if ( $IPfile =~ m/\.gz$/ ) {
    if ( $no_gunzip == 1 ) {
      print STDERR "[INFO] Unzipping is not possible.\n";
      exit 1;
    }
    open( $fh, "gunzip -c $IPfile |" );
  } else {
    open( $fh, $IPfile );
  }
  while ( my $line = <$fh> ) {
    chomp $line;
    next if ( $line =~ m/^#/ );
    my @F = split( /\t/, $line );
    $IP{"$F[0]:$F[5]:$F[2]"} = 1;
  }
}

my %clusters = ();
my %map      = ();
my @libsize  = ();
if ( !$pipe ) {
  foreach ( 0 .. $#files ) {
    push @libsize, 0;
  }
} else {

  # just one needed
  push @libsize, 0;
}

my $fh;
if ( $clusterfile =~ m/\.gz$/ ) {
  if ( $no_gunzip == 1 ) {
    print STDERR "[INFO] Unzipping is not possible.\n";
    exit 1;
  }
  open( $fh, "gunzip -c $clusterfile |" );
} else {
  open( $fh, $clusterfile );
}
while ( my $line = <$fh> ) {
  chomp $line;
  next if ( $line =~ m/^#/ );
  my @F = split( /\t/, $line );
  _checkBED( \@F, $clusterfile );
  if ( $F[4] !~ m/^\d+$/ ) {
    print STDERR "[ERROR] Wrong support level annotation.\n";
    print STDERR "[ERROR] Offending line: '$line'.\n";
    exit;
  }
  next if ( $F[4] < $minLevel );

  my @tmp = ();
  if ( !$pipe ) {
    foreach ( 0 .. $#files ) {
      push @tmp, 0;
    }
  } else {

    # just one needed
    push @tmp, 0;
  }
  $clusters{ $F[3] } = \@tmp;
  $map{"$F[0]:$F[5]"} = {} if ( not defined $map{"$F[0]:$F[5]"} );
  foreach my $pos ( $F[1] + 1 .. $F[2] ) {
    $map{"$F[0]:$F[5]"}->{$pos} = $F[3];
  }
}
close($fh);

if ( !$pipe ) {
  foreach my $idx ( 0 .. $#files ) {
    my $fh_sample;
    if ( $files[$idx] =~ m/\.gz$/ ) {
      open( $fh_sample, "gunzip -c $files[$idx] |" );
    } else {
      open( $fh_sample, $files[$idx] );
    }

    print STDERR "[INFO] Processing $files[$idx]\n";
    my $cnt_all      = 0;
    my $cnt_valid    = 0;
    my $cnt_assigned = 0;
    my $cnt_IP       = 0;
    if ( defined $savedRejected ) {
      open( OUT, ">$names[$idx].rejected.bed" );
    }
    while ( my $line = <$fh_sample> ) {
      my ( $type, $readcount ) =
        _checkRead( $line, $files[$idx], $scoreAsReadCount, $chrM, $IPfile, \%IP, \%map );
      $cnt_all += $readcount;
      if ( $type eq 'chrM' ) {

        # do not count to valid reads
      } elsif ( $type eq 'IP' ) {
        $cnt_valid += $readcount;
        $cnt_IP    += $readcount;
      } elsif ( $type eq 'nomatch' ) {
        $cnt_valid += $readcount;
        if ( defined $savedRejected ) {
          print OUT $line;
        }
      } else {
        $cnt_valid               += $readcount;
        $cnt_assigned            += $readcount;
        $clusters{$type}->[$idx] += $readcount;
      }
    }
    if ( defined $savedRejected ) {
      close(OUT);
    }
    $libsize[$idx] = $cnt_valid;
    my $pc_assigned = sprintf( "%.2f", $cnt_assigned / $cnt_valid * 100 );
    print STDERR "[INFO] A total of ", _formatNumber($cnt_all),   " reads screened.\n";
    print STDERR "[INFO] A total of ", _formatNumber($cnt_valid), " reads were considered.\n";
    print STDERR "[INFO] ",            _formatNumber($cnt_assigned),
      " ($pc_assigned%) reads were assigned to clusters.\n";
    if ($IPfile) {
      my $pc_IP = sprintf( "%.2f", $cnt_IP / $cnt_valid * 100 );
      print STDERR "[INFO] ", _formatNumber($cnt_IP),
        " ($pc_IP%) reads are classified as internal priming.\n";
    }
    close($fh_sample);
  }
} else {
  my $cnt_all      = 0;
  my $cnt_valid    = 0;
  my $cnt_assigned = 0;
  my $cnt_IP       = 0;
  while ( my $line = <STDIN> ) {
    my ( $type, $readcount ) =
      _checkRead( $line, 'STDIN', $scoreAsReadCount, $chrM, $IPfile, \%IP, \%map );
    $cnt_all += $readcount;
    if ( $type eq 'chrM' ) {

      # do not count to valid reads
    } elsif ( $type eq 'IP' ) {
      $cnt_valid += $readcount;
      $cnt_IP    += $readcount;
    } elsif ( $type eq 'nomatch' ) {
      $cnt_valid += $readcount;
      if ( defined $savedRejected ) {
        print OUT $line;
      }
    } else {
      $cnt_valid            += $readcount;
      $cnt_assigned         += $readcount;
      $clusters{$type}->[0] += $readcount;
    }
  }
  $libsize[0] = $cnt_valid;
  my $pc_assigned = sprintf( "%.2f", $cnt_assigned / $cnt_valid * 100 );
  print STDERR "[INFO] A total of ", _formatNumber($cnt_all),   " reads screened.\n";
  print STDERR "[INFO] A total of ", _formatNumber($cnt_valid), " reads were considered.\n";
  print STDERR "[INFO] ",            _formatNumber($cnt_assigned),
    " ($pc_assigned%) reads were assigned to clusters.\n";
  if ($IPfile) {
    my $pc_IP = sprintf( "%.2f", $cnt_IP / $cnt_valid * 100 );
    print STDERR "[INFO] ", _formatNumber($cnt_IP),
      " ($pc_IP%) reads are classified as internal priming.\n";
  }
}

if ( $#names == -1 ) {
  @names = ('clusterid');
  foreach my $idx ( 0 .. $#files ) {
    push @names, "sample" . ( $idx + 1 );
  }
} else {
  unshift @names, 'clusterid';
}

print join( "\t", @names ), "\n";
foreach my $cluster ( sort keys %clusters ) {
  my @tmp      = ();
  my $non_zero = 0;
  foreach my $idx ( 0 .. $#files ) {
    $non_zero++ if ( $clusters{$cluster}->[$idx] > 0 );
    if ($tpm) {
      push @tmp, sprintf( "%.5f", $clusters{$cluster}->[$idx] / $libsize[$idx] * 1000000 );
    } else {
      push @tmp, $clusters{$cluster}->[$idx];
    }
  }
  next if ( $non_zero == 0 );

  print "$cluster\t", join( "\t", @tmp ), "\n";
}

sub _checkBED {
  my @F    = @{ $_[0] };
  my $file = $_[1];
  if ( not defined $F[5] ) {
    print STDERR "[ERROR] Content of $file does not seem to be in BED format.\n";
    print STDERR "[ERROR] Offending line: '", join( "\t", @F ), "'.\n";
    exit;
  } else {
    if ( $F[5] !~ m/^(\+|-)$/ ) {
      print STDERR "[ERROR] Content of $file does not seem to be in BED format.\n";
      print STDERR "[ERROR] Offending line: '", join( "\t", @F ), "'.\n";
      exit;
    }
  }
}

sub _checkRead {
  my $lineL             = $_[0];
  my $fileL             = $_[1];
  my $scoreAsReadCountL = $_[2];
  my $chrML             = $_[3];
  my $IPfileL           = $_[4];
  my $IPL               = $_[5];
  my $mapL              = $_[6];

  chomp $lineL;
  my @F = split( /\t/, $lineL );
  _checkBED( \@F, $fileL );
  my $readCount = ($scoreAsReadCountL) ? $F[4] : 1;

  # determine the 3' end of the read
  my $pos;
  if ( abs( $F[2] - $F[1] ) == 1 ) {
    $pos = $F[2];
  } else {
    if ( $F[5] eq '+' ) {
      $pos = $F[2];
    } else {
      $pos = $F[1] + 1;
    }
  }

  # check if it maps to the mitochondrial genome
  if ($chrML) {
    return ( 'chrM', $readCount ) if ( $F[0] eq $chrML );
  }

  # check if it seem to come from IP
  if ($IPfileL) {
    if ( defined $IPL->{"$F[0]:$F[5]:$pos"} ) {
      return ( 'IP', $readCount );
    }
  }

  # check if it matches a cluster
  if ( defined $mapL->{"$F[0]:$F[5]"} ) {
    if ( defined $mapL->{"$F[0]:$F[5]"}->{$pos} ) {
      return ( $mapL->{"$F[0]:$F[5]"}->{$pos}, $readCount );
    }
  }

  return ( 'nomatch', $readCount );
}

sub _formatNumber {
  my $num = $_[0];
  $num = reverse $num;    # reverse the number's digits
  $num =~ s/(\d{3})/$1,/g;    # insert comma every 3 digits, from beginning
  $num = reverse $num;        # Reverse the result
  $num =~ s/^\,//;            # remove leading comma, if any

  return $num;
}
