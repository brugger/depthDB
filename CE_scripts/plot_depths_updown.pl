#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (27 Mar 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

my $norm = 10000000;

my @depths;

while ( my $file = shift ) {
  
  my $sample = $file;
  $sample =~ s/\..*//;

  my $i = 0;
  open ( my $in, $file) || die "Could not open '$file': $!\n";
  while (<$in>) {
    chomp;
    if ( $_ =~ /^===/) {
      $i=0;
      next;
    }
    my @F = split("\t");
    push @{$depths[$i++]}, $F[2];
  }


}

#print Dumper( \@depths );

my $i = 1;
#@depths = reverse @depths;
foreach my $pos ( @depths ) {
  print join("\t", $i++, mean(@$pos) ) . "\n";
}



# 
# 
# 
# Kim Brugger (27 Mar 2013)
sub mean {
  my (@values) = @_;

#  \''{s+=$1;s2+=($1*$1)} END {print sqrt((NR*s2-s*s)/(NR*(NR-1)))}'\'


  my $sum = 0;
  my $ss  = 0;
  map {$sum += $_; $ss += $_*$_} @values;

  my $NR = int(@values);

  my $mean = sprintf("%.2f", $sum/$NR);
  my $sd = sprintf("%.2f", sqrt(($NR*$ss-$sum*$sum)/($NR*($NR-1))));

  return ($mean);
  return ($mean-2*$sd, $mean, $mean+2*$sd);
}
