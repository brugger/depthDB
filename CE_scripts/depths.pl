#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Mar 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $bam = shift;
my $bed = shift;
$bed ||= "/data/refs/BRCA/BRCA12_exons.bed";


open (my $BED, $bed) || die "Could not open '$bed': $!\n";
while(<$BED>) {
  chomp;
  my @F = split("\t", $_);
  next if ( ! $F[2]);

  my $cmd = "/software/bin/samtools depth -r $F[0]:$F[1]-$F[2] $bam | ";
  open (my $p, $cmd);
  while(<$p>) {
    print;
  }
  close ($p);

  print "===\n";

}
