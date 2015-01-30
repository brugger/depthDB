#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (15 Oct 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $UKGTN    = shift;
my $TruSight = shift;

my %ukgtn;

open( my $in, $UKGTN);
while (<$in>) {
  chomp;
  $ukgtn{ $_ } ++;
}
close( $in);

open( $in, $TruSight );
while(<$in>) {
  
  next if (/#/);
  chomp;

  my @F = split("\t");
  $F[3] =~ s/^(.*?)\..*/$1/;
#  print "$_\n" if ( $ukgtn{ $F[3]});
  $ukgtn{ $F[3]}++;
}

foreach my $key (keys %ukgtn) {
  next if ( $ukgtn{ $key} > 1);
  
  print "$key\n";
}
  
