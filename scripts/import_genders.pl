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
open (my $in, $infile ) || die "Could not open '$infile': $!\n";
$infile =~ s/^\.\///;
$infile =~ s/\..*//;
my $cwd      = `pwd`;
chomp($cwd);


while (<$in> ) {
  chomp;
  my ( $sample_name, $mean_het, $mean_homo, $het_homo, $het_homo_score, $X_homo_het_ratio, $gender ) = split("\t", $_);
  my $sid = CTRU::depthDB::fetch_sample_id( $sample_name );
  print "Sample: $sample_name \n";
  next if ( ! $sid );

  $gender = 'F' if ($gender eq 'female');
  $gender = 'M' if ($gender eq 'male'  );

  print "$sid -- $gender \n";

  
  CTRU::depthDB::add_sample_stats( $sid, undef, undef, undef, undef, $gender);

}



