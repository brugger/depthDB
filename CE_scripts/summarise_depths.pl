#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (21 Jun 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;


my %min_depth  = ();
my %min_reads  = ();
my %exon_count = ();
my %genes      = ();

while ( my $file = shift ) {
  
#  print "$file\n";

  open( my $f, $file);
  while (<$f>) {

    chomp;
    my @F=split("\t");

    $F[2] ||= 0;
    
    my ($gene,$exon_nr) = split("_", $F[1]);
    $exon_nr =~ s/Exon//;

    $genes{ $gene }++;

    if ( ! $min_depth{ $file }{$gene} || 
	 $min_depth{ $file }{$gene} > $F[2]) {
      
      $min_depth{ $file }{$gene} = $F[2];
      $min_reads{ $file }{$gene} = $F[3];
    }
    if ( ! $exon_count{ $gene } || 
	 $exon_count{ $gene } < $exon_nr) {
      $exon_count{ $gene } = $exon_nr;
    }
    
  }
}


#print Dumper( \%min_depth);
#print Dumper( \%exon_count);

my @files = sort keys %min_depth;

my @samples = @files;
map { s/\..*//} @samples;

print join("\t", "Gene", @samples) . "\n";

foreach my $gene ( sort keys %genes ) { 
  my @values;
  push @values, $gene;
  push @values, $exon_count{ $gene } ;

  foreach my $file ( @files  ) {
    push @values, $min_depth{ $file}{ $gene }, $min_reads{ $file }{ $gene } ;
  }
  
  print join("\t", @values ) . "\n";
    
}
