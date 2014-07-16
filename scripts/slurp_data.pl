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
$infile =~ s/\..*//;
my $cwd      = `pwd`;
chomp($cwd);



my $sample_sequence_name = $infile;
$sample_sequence_name =~ s/.*\///;
my $sample_name = substr($sample_sequence_name, 0, 7);

my $plate_name = 'Unknown';
$plate_name = $1 if ($infile =~ /(CP\d+)/ || $cwd =~ /(CP\d+)/);
print " $sample_name => $sample_sequence_name, $plate_name\n";

my $sid = CTRU::depthDB::add_sample( $sample_name );

readin_csv("$infile.var.csv");
readin_stats("$infile.bam.flagstat", "$infile.bam.isize");


# 
# 
# 
# Kim Brugger (04 Dec 2013)
sub readin_stats {
  my ($flagstat_file, $isize_file) = @_;

  my %res;
  open (my $in, $flagstat_file) || die "Could not open '$flagstat_file': $!\n";
  while(<$in>) {
    print;
    if ( /(\d+) .*total/) {
      $res{total_reads} = $1;
    }
    elsif ( /^(\d+) .*duplicates/) {
      $res{dup_reads} = $1;
    }
    elsif ( /(\d+) .*mapped \((.*?)%/) {
      $res{mapped_reads} = $1;
      $res{mapped_perc} = $2;
    }
    elsif ( /(\d+) .*properly paired/) {
      $res{properly_paired} = $1;
    }
    elsif ( /(\d+) .*singletons/) {
      $res{singletons} = $1;
    }
  }
  close( $in );

  open(  $in, $isize_file) || die "Could not open '$isize_file': $!\n";
  while(<$in>) {
    chomp;
    if ( /^MEDIAN_INSERT_SIZE/) {
      my @rows = split("\t");
      $_ = <$in>;
      chomp;
      my @fields = split("\t");

      for( my $i=0;$i< @rows;$i++) {
	$res{ lc($rows[ $i ])} = $fields[ $i ];
      }
      last;
    }
  }
  close( $in );


#  print Dumper(\%res);

  
  CTRU::depthDB::add_sample_sequence_stats( $ssid, $res{total_reads}, $res{mapped_reads}, $res{dup_reads}, $res{mean_insert_size});

}



