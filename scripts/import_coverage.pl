#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (04 Dec 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib '/software/packages/depthDB/modules';
use CTRU::depthDB;

my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';
my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");

use Getopt::Std;

my $infile = shift || usage();


my $sample_name = $infile;
$sample_name =~ s/.*\///;
$sample_name =~ s/\_mem.*//;


my $sid = CTRU::depthDB::fetch_sample_id( $sample_name );

print " $infile :: $sample_name $sid\n";

my @coverages;
open(my $in, $infile ) || die "Could not open '$infile': $!\n";
while (<$in>) {
  chomp;

  my ($rid, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19) = split("\t", $_);


#  print "[$sid] -- $rid, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19 \n";
  push @coverages, [$sid, $rid, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19];

#  last;
}
  
  
CTRU::depthDB::add_coverages( @coverages );
