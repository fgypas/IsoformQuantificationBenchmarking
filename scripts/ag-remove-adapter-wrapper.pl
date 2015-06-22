use strict;
use warnings;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;

my $aln2seq_IOeff;
my $blastscores;
my $removeAdaptor;
my $adapter;
my $help;
my $debug;
my $bin = 5000000;

my $result = GetOptions(
  "aln2seq=s"       => \$aln2seq_IOeff,
  "blastscores=s"   => \$blastscores,
  "removeAdaptor=s" => \$removeAdaptor,
  "adapter=s"       => \$adapter,
  "a=s"             => \$adapter,
  "binsize=i"       => \$bin,
  "b=i"             => \$bin,
  "help"            => \$help,
  "debug"           => \$debug,
  "h"               => \$help
);

my $showhelp = 0;
$showhelp = 1 if ( not defined $ARGV[0] );
$showhelp = 1 if ( not defined $adapter );
$showhelp = 1 if ( not defined $aln2seq_IOeff );
$showhelp = 1 if ( not defined $blastscores );
$showhelp = 1 if ( not defined $removeAdaptor );

if ($showhelp) {
  print STDERR "Usage $0 --aln2seq=path --blastscores=path --removeAdaptor=path ";
  print STDERR "-a ACGTGACGA -b 5000000 file.fa.gz\n\n";
  exit;
}

if ( not -e $ARGV[0] ) {
  print STDERR "[ERROR] File '$ARGV[0]' not found.\n";
  exit;
}

if ( $ARGV[0] !~ m/fa.gz$/ ) {
  print STDERR "[ERROR] File '$ARGV[0]' does not look like gzipped FASTA file.\n";
  exit;
}

if ( not -e $aln2seq_IOeff ) {
  print STDERR "[ERROR] Please check path to aln2seq_IOeff executable.\n";
  print STDERR "[ERROR] $aln2seq_IOeff\n";
  exit;
}

if ( not -e $blastscores ) {
  print STDERR "[ERROR] Please check path to blast.scores file.\n";
  exit;
}

# make a temp file
my $outdir = '.';
if ( $ARGV[0] =~ m/(.*\/)(.*)/ ) {
  $outdir = $1;
}

# write adapter to disk
my @files = ( File::Temp->new( DIR => $outdir, SUFFIX => '.adapter' ) );
$files[0]->print("$adapter\n");
$files[0]->eof();

# convert gzipped fasta file to sol file
push @files, File::Temp->new( DIR => $outdir, SUFFIX => '.sol' );
my $cnt   = 0;
my $tmpfh = $files[$#files];
open( IN, "gunzip -c $ARGV[0] |" ) || die "can't open pipe to $ARGV[0]";
while (<IN>) {
  chomp;
  next if ( $_ =~ m/^>/ );
  $cnt++;
  $tmpfh->print("$_\t1\n");
  if ( $cnt % $bin == 0 ) {
    $tmpfh->eof();
    push @files, File::Temp->new( DIR => $outdir, SUFFIX => ".sol" );
    $tmpfh = $files[$#files];
  }
}
close(IN);
$tmpfh->eof();

foreach my $file (@files) {
  my $in = $file->filename;
  next if ( $in !~ m/.sol$/ );
  my $out = File::Temp->new( DIR => $outdir, SUFFIX => ".sol.aln" );
  push @files, $out;
  my $outfn     = $out->filename;
  my $adapterfn = $files[0]->filename;
  `$aln2seq_IOeff $in $adapterfn $blastscores > $outfn`;
  if ($debug) {
    open( D, $outfn );
    while (<D>) {
      print STDERR "$_\n";
    }
    close(D);
  }
}

# remove adaptor
$cnt = 0;
foreach my $file (@files) {
  my $in = $file->filename;
  next if ( $in !~ m/.sol.aln$/ );
  open( IN2, "$removeAdaptor $in |" ) || die "can't open pipe";
  while (<IN2>) {
    chomp;
    my @F = split(/\t/);
    next if ( length( $F[0] ) < 15 );
    $cnt++;
    print ">$cnt\n$F[0]\n";
  }
  close(IN2);
}
