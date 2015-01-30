#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (28 Apr 2010), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

use lib '/software/packages/ensembl/modules/';
use lib '/software/packages/bioperl/';

# Sets up dynamic paths for EASIH modules...
# Makes it possible to work with multiple checkouts without setting 
# perllib/perl5lib in the enviroment.
BEGIN {
  my $DYNAMIC_LIB_PATHS = 0;
  if ( $DYNAMIC_LIB_PATHS ) {
    my $path = $0;
    if ($path =~ /.*\//) {
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules");
      $path =~ s/(.*)\/.*/$1/;
      push @INC, "$path/modules" if ( -e "$path/modules" && ! grep /^$path\/modules/, @INC);
    }
    else {
      push @INC, "../modules" if ( -e "../modules");
      push @INC, "./modules" if ( -e "./modules");
    }
  }
  else {
    push @INC, '/home/kb468/easih-toolbox/modules/';
  }
}


use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;

my $species     = "human";
my $host        = 'mgsrv01.medschl.cam.ac.uk';
my $user        = "easih_ro";
# get registry
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(-host => $host,-user => $user);

my $exit_count = 9;
# get adaptors
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');
my $ta  = $reg->get_adaptor($species, 'core', 'transcript');

my $transcript_id  = 'NM_080680';

while ($transcript_id = <>) {
  
  chomp ($transcript_id );
#  print "$transcript_id";

  my @genes = @{ $ga->fetch_all_by_external_name( $transcript_id )};

  my %exons;

  my $gene_name = "";
  my $refseq_id = "";
  my $ccds_id   = "";
  foreach my $gene ( @genes ) {
    

    my $seq_region = undef;
    my $strand    = 1;
    my $transcripts = $gene->get_all_Transcripts();
    while ( my $transcript = shift @{$transcripts} ) {
      my @db_entries = @{ $transcript->get_all_DBEntries('RefSeq_dna')};

      my @dbentries = @{ $transcript->get_all_DBEntries() };
      
      $gene_name = $gene->external_name();
      next if ( $gene_name ne $transcript_id);
      
#      print "----------------\n";

      my ($refseq_entry, $ccds_entry) = (0,0);
      
      foreach my $db_ent ( @dbentries ) {
#	print join("\t", $db_ent->display_id() , $db_ent->dbname() ) ."\n";
	$refseq_entry++ if ( $db_ent->dbname() eq "RefSeq_dna");
	$ccds_entry++   if ( $db_ent->dbname() eq "CCDS");
	$ccds_id = $db_ent->display_id() if ( $db_ent->dbname() eq "CCDS");
      }
      
      next if ( ! $refseq_entry || ! $ccds_entry);

      foreach my $db_entry (@db_entries) {
	
	$refseq_id = $db_entry->primary_id();
	foreach my $exon ( @{ $transcript->get_all_translateable_Exons() } ) {
	  
	  $seq_region = $exon->slice->seq_region_name();
	  my $start      = $exon->start();
	  my $end        = $exon->end();
	  ($start, $end) = ($end, $start) if ( $start > $end );
	  $strand     = $exon->strand();
	  
	  $exons{$start} = $end;
	}
      }
    }

    next if ( ! $seq_region );
    
    
    my $counter = 1;
    $counter = keys %exons  if ( $strand == -1 );
    foreach my $start ( sort  { $a <=> $b } keys %exons) {
      print join("\t", $seq_region, $start, $exons{$start} , $gene_name ."_Exon". $counter) . "\n";
      $counter++ if ( $strand == 1 );
      $counter-- if ( $strand == -1 );
    }
  }
}

# 
# 
# 
# Kim Brugger (09 Nov 2010)
sub usage {
  
  $0 =~ s/.*\///;

#  print "USAGE: $0 -b[am file] -i[indel vcf file] -s[np vcf file] -T<ranform, use if mapped against hg18> -B[ait file] -l[eeway, default 100 bp] -c[ount bases, need a -b as well]\n";
  print "USAGE: $0 -b[am file] -v[ariant vcf file] -T<ranform, use if mapped against hg18> -B[ait file] -l[eeway, default 100 bp] -c[ount bases, need a -b as well]\n";

#  print "\nor extrapolate the standard <bam, SNP vcf, indel vcf, output files> with the -Q <basefile name> option\n";
  print "\nor extrapolate the standard <bam, vcf, output files> with the -Q <basefile name> option\n";
  print "EXAMPLE: $0 -Q [base name] -T<ransform>\n";
  print "\n";

  
  print "USAGE: -o[output file]\n";  

  exit;
}



