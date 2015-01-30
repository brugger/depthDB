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
$infile =~ s/^\.\///;
$infile =~ s/\..*//;
my $cwd      = `pwd`;
chomp($cwd);


my $qc_file = "$infile.vcf.qc";
if ( ! -e "$qc_file" ) {
  $qc_file =~ s/.*\///;
  if ( -e "vcfs/$qc_file" ) {
    $qc_file = "vcfs/$qc_file";
  }
  else { 
    die "Cannot find qc file (vcfs/$qc_file)!\n";
  }
}


my $sample_name = $infile;
$sample_name =~ s/.*\///;
$sample_name =~ s/.bam//;
$sample_name =~ s/\_mem//;

print " $infile :: $sample_name\n";

my $sid = CTRU::depthDB::add_sample( $sample_name );

readin_stats("$infile.bam.flagstat", "$infile.bam.isize", $qc_file);





# 
# 
# 
# Kim Brugger (04 Dec 2013)
sub readin_stats {
  my ($flagstat_file, $isize_file, $vcf_qc_file) = @_;

  my %res;
  open (my $in, $flagstat_file) || die "Could not open '$flagstat_file': $!\n";
  while(<$in>) {
#    print;
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



  my $gender = 'U';

  open(  $in, $vcf_qc_file) || die "Could not open '$vcf_qc_file': $!\n";
  while (<$in> ) {
#    print ;
    chomp;
    next if ( /^Sample/);
    my @F = split("\t", $_ );
    
#    my ( $sample_name, $mean_het, $mean_homo, $het_homo, $het_homo_score, $X_homo_het_ratio, $gender ) = split("\t", $_);

#    print "G: $gender $F[0] \n";

    $gender = 'F' if ($F[6] eq 'female');
    $gender = 'M' if ($F[6] eq 'male'  );
#    print "G: $gender $F[0] \n";
    last;
  }
  close( $in );

  
  CTRU::depthDB::add_sample_stats( $sid, $res{total_reads}, $res{mapped_reads}, $res{dup_reads}, $res{mean_insert_size}, $gender);

}



