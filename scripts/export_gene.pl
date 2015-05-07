#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (07 May 2015), contact: kbr@brugger.dk

use strict;
use warnings;
use Data::Dumper;


use lib '/software/packages/depthDB/modules';
use CTRU::depthDB;

my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';
my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");


my %bed_regions;

while( my $line = <>) {
  chomp $line;

  my($gene, $transcript) = split(/\s/, $line);


  my $regions = CTRU::depthDB::fetch_regions_by_gene( $gene );
  

  foreach my $region ( @$regions ) {

    next if ( $transcript && $$region{ 'refseq' } !~ /$transcript/);
#    print Dumper( $region );

    my $region_hash = CTRU::depthDB::fetch_region_hash( $$region{ 'rid' });

    $bed_regions{ $$region_hash{ 'chr'} }{ $$region_hash{ 'start'}}{ $$region_hash{ 'end'}}++;

  }
}
  
  
foreach my $chr ( sort keys %bed_regions ) {
  foreach my $start ( sort keys %{$bed_regions{ $chr }}) {
    foreach my $end ( sort keys %{$bed_regions{ $chr }{$start}}) {
      
      print join("\t", $chr, $start, $end )."\n";
    }
  }
}
