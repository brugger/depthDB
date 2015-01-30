#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (19 Dec 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my %res;

while(<>)  {
  chomp;
  my ($gene, $refseq, $ccds, $pos, $name, $low, $reads) = split("\t");
  $res{ $gene } { $refseq} { $ccds } {exons}++;
  $res{ $gene } { $refseq} { $ccds } {exons_fail}++ if ( $low < 20);
  

  if ( !$res{ $gene } { $refseq} { $ccds } { low }) {
    $res{ $gene } { $refseq} { $ccds } { low }   = $low;
    $res{ $gene } { $refseq} { $ccds } { reads } = $reads;
  }

  if ( $res{ $gene } { $refseq} { $ccds } { low } && $res{ $gene } { $refseq} { $ccds } { low } > $low ) {
    $res{ $gene } { $refseq} { $ccds } { low }   = $low;
    $res{ $gene } { $refseq} { $ccds } { reads } = $reads;
  }
}


#print Dumper( \%res );

foreach my $gene ( sort keys %res ) {
  foreach my $refseq ( sort keys %{$res{ $gene }} ) {
    foreach my $ccds ( sort keys %{$res{ $gene }{ $refseq }} ) {
      print join("\t", $gene, $refseq, $ccds, 
		 $res{ $gene } { $refseq} { $ccds } { low },
		 $res{ $gene } { $refseq} { $ccds } { reads }) . "\n";
    }
  }
}


