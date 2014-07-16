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

my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';
my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");


use Getopt::Std;
my $opts = 's:';
my %opts;
getopts($opts, \%opts);

my $infile = $opts{s};



my $sample_sequence_name = $infile;
$sample_sequence_name =~ s/.*\///;
$sample_sequence_name =~ s/\..*//;
my $plate_name = 'CP0054';

my $sample_name = substr($sample_sequence_name, 0, 7);

print " $sample_name => $sample_sequence_name\n";


open( my $in, $infile) || die "Could not open '$infile': $!\n";
while(<$in>) {
  next if (!/#  capture:/);
  chomp;
  s/#  capture: //;
  my ($name, $min, $max, $mean, $lows, $missing, $transcript) = split("\t");

  print "  my ($name, $min, $max, $mean, $lows, $missing, $transcript) = \n";

  my $rid = CTRU::depthDB::fetch_region_id_by_name( $name );

  if ( !$rid ) {
    print STDERR "could not find region: '$name', skipping on to the next one\n";
    next;
  }

  my $cid = CTRU::depthDB::add_coverage($ssid, $rid, $min, $mean, $max, $lows, $missing);
  
  
  
  
}
