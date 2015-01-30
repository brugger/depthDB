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
use EASIH::DB;

my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';
my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");

use Getopt::Std;


my $sample_name = shift;
my $sid = CTRU::depthDB::fetch_sample_id( $sample_name );

if (!$sid) {
  print "$sample_name not present in the database\n";
  exit -1;
}

print "$sample_name -- SID :: $sid\n";
print "==============================================\n";
print "                 WARNING                     \n";
print "This will purge $sample_name from the database\nYou have 5 sec to press ctrl-c to abort\n";
print "==============================================\n";


sleep (5);

my $delete_coverage = "delete from coverage where sid = '$sid';";
my $delete_sample   = "delete from sample where sid = '$sid';";

EASIH::DB::do( $dbi, $delete_coverage );
EASIH::DB::do( $dbi, $delete_sample   );


