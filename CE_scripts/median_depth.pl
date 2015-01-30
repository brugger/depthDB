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

my @all_depths;

open (my $BED, $bed) || die "Could not open '$bed': $!\n";
while(<$BED>) {
  chomp;
  my @F = split("\t", $_);
  next if ( ! $F[2]);

  $F[1] -= 5;
  $F[2] += 5; 
  
  median($F[3], "/software/bin/samtools depth -Q 20 -r $F[0]:$F[1]-$F[2] $bam | ");   

}

@all_depths = sort {$a <=> $b} @all_depths;

print "Overall\t" . $all_depths[ int(@all_depths/2) ] ."\n";

# 
# 
# 
# Kim Brugger (22 Mar 2013)
sub median {
  my ($id, $cmd) = @_;

  my @depths;

  open (my $p, $cmd);
  while(<$p>) {
    chomp;
    my @D = split("\t", $_);
    push @all_depths, $D[2];
    push @depths, $D[2];
  }

  push @depths, 0 if (@depths == 0);


  @depths = sort {$a <=> $b} @depths;
  
  print join("\t", "$id", $depths[0],  $depths[ int(@depths/2) ], $depths[-1]) ."\n";

  
}

