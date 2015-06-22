use strict;
use warnings;
use Getopt::Long;

# TODO:
# option for output file(s)
# if possible, produce only 1 output file including both unique and multi mappers

my $newfile1;
my $newfile2;
my $samToBam;
my $samToBed;
my $help;
my $transcriptome;
my $bedToGenome;
my $mapfile;
my $if1file;
my $threads = 1;
my @exclude = ();

GetOptions(
  "bed-unique:s"  => \$newfile1,
  "bed-multi:s"   => \$newfile2,
  "samtools:s"    => \$samToBam,
  "samToBed:s"    => \$samToBed,
  "bedToGenome:s" => \$bedToGenome,
  "transcriptome" => \$transcriptome,
  "mapfile:s"     => \$mapfile,
  "if1file:s"     => \$if1file,
  "threads:i"     => \$threads,
  "exclude:s"     => \@exclude,
  "help"          => \$help,
  "h"             => \$help
);

my $showhelp = 0;
$showhelp = 1 if ( defined $help );
$showhelp = 1 if ( not defined $samToBam );
$showhelp = 1 if ( not defined $samToBed );
$showhelp = 1 if ( not defined $ARGV[0] );
if ( defined $transcriptome ) {
  $showhelp = 1 if ( not defined $mapfile );
  $showhelp = 1 if ( not defined $if1file );
}

if ($showhelp) {
  print STDERR "Usage: $0 --samtools=path --bamToBed=path --bed-unique=path file.sam\n";
  print STDERR "\n";
  exit;
}

if ( not defined $newfile1 ) {
  print STDERR "Usage: $0 --samtools=path --bamToBed=path --bed-unique=path file.sam\n";
  print STDERR "[ERROR] Required output filepath ('--bed-unique') missing.\n";
  exit 1;
}

if ( not -e $ARGV[0] ) {
  ( my $tmpbam = $ARGV[0] ) =~ s/\.sam$/.bam/;
  if ( -e $tmpbam ) {
    print STDERR "[INFO] Found BAM file. Skipping conversion.\n";
  } else {
    print STDERR "[ERROR] File '$ARGV[0]' not found.\n";
    exit 1;
  }
}
if ( not -e $samToBam ) {
  print STDERR "[ERROR] File '$samToBam' not found.\n";
  exit 1;
}
if ( not -e $samToBed ) {
  print STDERR "[ERROR] File '$samToBed' not found.\n";
  exit 1;
}

if ($transcriptome) {
  if ( not -e $bedToGenome ) {
    print STDERR "[ERROR] File '$bedToGenome' not found.\n";
    exit 1;
  }
  if ( not -e $mapfile ) {
    print STDERR "[ERROR] File '$mapfile' not found.\n";
    exit 1;
  }
  if ( not -e $if1file ) {
    print STDERR "[ERROR] File '$if1file' not found.\n";
    exit 1;
  }
} else {
  if ( not defined $newfile2 ) {
    print STDERR "Usage: $0 --samtools=path --bamToBed=path --bed-unique=path --bed-multi=path file.sam\n";
    print STDERR "[ERROR] Required output filepath ('--bed-multi') missing.\n";
    exit 1;
  }
}


my %excludeHash = ();
foreach my $e (@exclude) {
  $excludeHash{$e} = 1;
}

( my $bam = $ARGV[0] ) =~ s/\.sam$/.bam/;
if ( !-e $ARGV[0] and -e $bam ) {
  # do nothing
  print STDERR "[INFO] BAM file found. Not doing conversion.\n";
} elsif ( -e $ARGV[0] ) {
  print STDERR "[INFO] Converting SAM -> BAM.\n";
  `$samToBam view -S -b $ARGV[0] > $bam`;
} else {
  print STDERR "[ERROR] No SAM/BAM file.\n";
  exit;
}

print STDERR "[INFO] Converting BAM -> BED.\n";
( my $bed = $ARGV[0] ) =~ s/\.sam$/.bed/;

`$samToBam view $bam | python $samToBed -p $threads  > $bed`;

# for transcriptome alignments
if ($transcriptome) {
  # exclude reads that map to the minus strand of the transcriptome (this doesn't exist)
  my %tmphash = ();
  open( F, $bed );
  while (<F>) {
    chomp;
    my @F = split(/\t/);
    $tmphash{ $F[3] } = 1 if ( $F[5] eq '-' );
  }
  close(F);

  ( my $newbed = $bed ) =~ s/\.bed$/.genomic.bed/;
  open( F, $bed );
  open( O,
    "| python $bedToGenome -p $threads -i $if1file -d $mapfile > $newbed"
  );
  while (<F>) {
    chomp;
    my @F = split(/\t/);
    next if ( defined $tmphash{ $F[3] } );
    print O join( "\t", @F ), "\n";
  }
  close(O);
  close(F);

  `rm $bed`;
  $bed = $newbed;
  print STDERR "[INFO] $bed\n";

  my %h = ();
  open( F, $bed );
  while (<F>) {
    chomp;
    my @F   = split(/\t/);
    my $key = "$F[0]:$F[1]:$F[2]:$F[5]:$F[6]:$F[7]";
    $h{ $F[3] } = {} if ( not defined $h{ $F[3] } );
    $h{ $F[3] }->{$key} = $F[4];
  }
  close(F);

  ( $newfile1 = $bed ) =~ s/\.bed$/.unique.bed/;
  open( O1, ">$newfile1" );
  foreach my $read ( keys %h ) {
    my @K = keys %{ $h{$read} };
    # consider only unique mappers
    if ( $#K == 0 ) {
      my @tmp = split( /:/, $K[0] );
      # if it matches over a splice junction
      if ( $tmp[5] =~ m/\|/ ) {
        my @N = split( /\|/, $tmp[5] );
        my $flag = 1;
	# if less than 5 nt on each side we do not consider it
        foreach my $n (@N) {
          $flag = 0 if ( $n < 5 );
        }
        $flag = 0 if ( defined $excludeHash{ $tmp[0] } );
        if ( $flag == 1 ) {
          print O1 "$tmp[0]\t$tmp[1]\t$tmp[2]\t$read\t$h{$read}->{$K[0]}\t$tmp[3]\t$tmp[4]\n";
        }
      }
    }
  }
  close(O1);

# for genome alignments
} else {
  my %h = ();
  open( F, $bed );
  while (<F>) {
    chomp;
    my @F = split(/\t/);
    next if ( defined $excludeHash{ $F[0] } );
    $h{ $F[3] }++;
  }
  close(F);

  print STDERR "[INFO] Splitting unique/multi mappers.\n";
  open( F,  $bed );
  open( O1, ">$newfile1" );
  open( O2, ">$newfile2" );
  while (<F>) {
    chomp;
    my @F = split(/\t/);
    next if ( defined $excludeHash{ $F[0] } );
    if ( $h{ $F[3] } == 1 ) {
      print O1 "$_\n";
    } else {
      print O2 "$_\n";
    }
  }
  close(O1);
  close(O2);
  close(F);
}

`rm $bed`;
if ( -e "$newfile1.gz" ) {
  `rm $newfile1.gz`;
}
if ( -e "$newfile2.gz" ) {
  `rm $newfile2.gz`;
}
`gzip $newfile1`;
if ( -e $newfile2 ) {
  `gzip $newfile2`;
}

if ( -e $ARGV[0] ) {
  `rm $ARGV[0]`;
}

