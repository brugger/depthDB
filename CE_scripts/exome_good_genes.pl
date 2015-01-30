#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (15 Oct 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my %genes;

while(<>) {
  chomp;
  my ( $exon, $low, $mean, $high) = split("\t");

  (my $gene, $exon) = split("_", $exon);

  $gene =~ s/\..*//;

  if ( ! $genes{ $gene }){
    $genes{ $gene } = $low;
  }
  elsif ( ! $genes{ $gene } > $low ){
    $genes{ $gene } = $low;
  }


}

my $MIN_DEPTH = 20;

foreach my $gene ( sort keys %genes) {
  print "$gene\t$genes{$gene}\n" if ($genes{$gene} < $MIN_DEPTH);
}
  
