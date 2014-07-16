#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (19 Nov 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/software/packages/depthDB/modules';
use CTRU::depthDB;

use Getopt::Std;
my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';

my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");
my $file = shift;
open( my $in, $file) || die "Could not open '$file': $!\n";

while(<$in>) {
  next if (/#/);
  next if (/^\z/);
  print;
  chomp;
  my ($chr, $start, $end, $name, $ref, $CCDS) = split("\t", $_);

  my $rid = CTRU::depthDB::fetch_region_id_by_position($chr, $start, $end);
  if ( !$rid ) {
    $rid = CTRU::depthDB::add_region($chr, $start, $end);
  }

  CTRU::depthDB::add_transcript($rid, $name, $ref, $CCDS) if ( $rid );
  print "RID :: $rid\n";
  
}
