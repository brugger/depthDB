#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (22 Mar 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my $bed = shift;

my $MIN_DEPTH = 30;

my @regions;

open (my $BED, $bed) || die "Could not open '$bed': $!\n";
while(<$BED>) {
  chomp;
  my @F = split("\t", $_);
  next if ( ! $F[2]);

  $F[1] -= 5;
  $F[2] += 5; 

  push @regions, \@F;

}

my %gene_stats;


while (my $bam = shift) {


  foreach my $region ( @regions ) {
    my ( $chr, $start, $end, $exon) = @{$region};

  
    my ($min, $median, $max) = min_median_max("/software/bin/samtools depth -r $chr:$start-$end $bam | ");   
    $min ||= 0;

#    print join("\t", $exon, $min, $median, $max) . "\n";
    
    my $gene = $exon;
    $gene =~ s/_.*//;

    $gene_stats{ $bam } { $gene } = $min if ( !$gene_stats{ $bam }{ $gene } || $gene_stats{ $bam }{ $gene } > $min);
#    $gene_stats{ $bam } { $gene } = $min if ( !$gene_stats{ $bam }{ $gene } || $gene_stats{ $bam }{ $gene } > $min);

  }
}

#print Dumper( \%gene_stats );

my $done_header = 0;

foreach my $bam ( sort keys %gene_stats ) {

  if ( ! $done_header ) {
    print join("\t", "Sample", sort keys %{$gene_stats{$bam}}) ."\n";
    $done_header++;
  }

  my @line;
  push @line, $bam;

  foreach my $gene ( sort keys %{$gene_stats{$bam}} ) {
    push @line, $gene_stats{$bam}{$gene};
  }

  print join("\t", @line)."\n";
}

# 
# 
# 
# Kim Brugger (22 Mar 2013)
sub min_median_max {
  my ($cmd) = @_;

  my @depths;

  open (my $p, $cmd);
  while(<$p>) {
    chomp;
    my @D = split("\t", $_);
    push @depths, $D[2];
  }

  @depths = sort {$a <=> $b} @depths;
  
  return ($depths[0],  $depths[ int(@depths/2) ], $depths[-1]);
  
}

