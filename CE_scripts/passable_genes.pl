#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (09 Oct 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my %gene_stats;

my $MIN_DEPTH = 20;

while( my $infile = shift) {

  open(my $in, $infile);
  while(<$in>) { 
    chomp;
    
    my ( $gene, $min, $median, $max) = split("\t");
    $min ||= 0;
    if ( $min < $MIN_DEPTH) {
      $gene_stats{ $gene }{$infile} = $min;
    }
  }
}

#print Dumper( \%gene_stats );
    

foreach my $exon ( sort keys %gene_stats ) {
  my @values;
  foreach my $file ( keys %{$gene_stats{ $exon}} ) {
    push @values, "$file:$gene_stats{ $exon}{$file}";
  }

  print join("\t", $exon, @values) ."\n";
  
}
