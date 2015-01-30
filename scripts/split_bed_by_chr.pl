#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (01 Sep 2014), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my %fhs;
my $infile = shift;
open( my $in, $infile ) || die "Could not open '$infile': $!\n";
while (<$in>) {
  chomp;
  my @F = split("\t");

  if (! $fhs{ $F[0] } ){
    open(my $fh, ">$F[0].bed") || die "Could not open '$F[0].bed': $!\n";
    $fhs{ $F[0] } = $fh;
  }
  my $fh = $fhs{ $F[0] };

  print $fh "$_\n";
}
  
