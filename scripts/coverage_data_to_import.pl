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
#$dbname = "frags_exome";
my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");

use Getopt::Std;

while (my $infile = shift) {
  chomp( $infile );

  my $sample = $infile;
  $sample =~ s/.bam.*//;
  $sample =~ s/_mem//;
  $sample =~ s/b(eta\d)\//b$1_/i;

  $sample =~ s/.*\///;
  $sample =~ s/\..*//;


  my $sid = CTRU::depthDB::fetch_sample_id( $sample );

  if ( ! $sid ) {
    die  "Uknown sample: '$sample' | $sid\n";
    next;
  }

#  print "$infile --> $sample --> $sid\n";

  open ( my $in, $infile ) || die "Could not open '$infile': $!";
  while( <$in> ) {
    next if (/No more bed/);
    print "$sid\t$_";
  }
  close( $in )
    

}

  
