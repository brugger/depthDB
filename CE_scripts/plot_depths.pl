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

my %reads = ('A1_S1'   => 121704352,
	     'B2_S2'   => 155137843,
	     'C3_S3'   => 35102065,
	     'D4_S4'   => 66119342,
	     'E5_S5'   => 90827124,
	     'F6_S6'   => 111287797,
	     'G7_S7'   => 41867993,
	     'H8_S8'   => 42986738,
	     'I9_S9'   => 85233235,
	     'J10_S10' => 20386913,
	     'K11_S11' => 117346512,
	     'L12_S12' => 69577773);

my @depths;

while ( my $file = shift ) {
  
  my $sample = $file;
  $sample =~ s/\..*//;

  my $i = 0;
  open ( my $in, $file) || die "Could not open '$file': $!\n";
  while (<$in>) {
    chomp;
    my @F = split("\t");
    push @{$depths[$i++]}, sprintf("%.2f", $F[2]/ $reads{ $sample }*$norm);
#    last if ( $i > 20 );
  }
  
}

#print Dumper( \@depths );

my $i = 1;
foreach my $pos ( @depths ) {
  print join("\t", $i++, mean(@$pos), @$pos ) . "\n";
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

  return ($mean-2*$sd, $mean, $mean+2*$sd);
}
