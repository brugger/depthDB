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

#use lib '/software/packages/ensembl_74/modules/';
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
#$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

my $exit_count = 9;
# get adaptors
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');
my $ta  = $reg->get_adaptor($species, 'core', 'transcript');


while (my $transcript_id = <>) {
  
  chomp ($transcript_id );
#  print "$transcript_id";

  my @genes = @{ $ga->fetch_all_by_external_name( $transcript_id )};
#  my @genes = @{ $ga->fetch_all_by_display_label( $gene_id )};

  if ( ! @genes ) {

    print "'$transcript_id' not found in the Ensembl database\n";
  }


#  print Dumper( @genes );

#  next ;

  my %exons;
  my %transcript_lengths;
  my %refseq2ccds;
  my $gene_name = "";
  foreach my $gene ( @genes ) {
    
    my $seq_region = undef;
    my $strand    = 1;
    my $transcripts = $gene->get_all_Transcripts();

    while ( my $transcript = shift @{$transcripts} ) {

      my @dbentries = @{ $transcript->get_all_DBEntries() };
      
      $gene_name = $gene->external_name();
#      print " $gene_name -- id\n";
#      next if ( $gene_name ne $gene_id);

      my ($refseq_entry, $ccds_entry) = (0,0);
      my ($refseq_id, $ccds_id)   = ("", "");
      
      foreach my $db_ent ( @dbentries ) {
#	print join("\t", $db_ent->display_id() , $db_ent->dbname() ) ."\n";
	$refseq_entry++ if ( $db_ent->dbname() eq "RefSeq_dna");
	$refseq_entry++ if ( $db_ent->dbname() eq "RefSeq_mRNA");
	$refseq_id = $db_ent->display_id() if ( $db_ent->dbname() eq "RefSeq_dna");
	$refseq_id = $db_ent->display_id() if ( $db_ent->dbname() eq "RefSeq_mRNA");
	$ccds_entry++   if ( $db_ent->dbname() eq "CCDS");
	$ccds_id = $db_ent->display_id() if ( $db_ent->dbname() eq "CCDS");
      }



      
#      print "$transcript_id -- $refseq_entry || ! $ccds_entry\n";
      next if ( ! $refseq_entry || ! $ccds_entry);
#      print " $refseq_id ne $transcript_id \n";

      next if ( $refseq_id ne $transcript_id);


      $refseq2ccds{ $refseq_id } = $ccds_id;
      my @db_entries = (@{$transcript->get_all_DBEntries('RefSeq_dna')}, @{ $transcript->get_all_DBEntries('RefSeq_mRNA')});


      foreach my $db_entry (@db_entries) {
	
#	print "$refseq_id == " . $db_entry->display_id() . "\n";

	next if ( $db_entry->display_id() ne $refseq_id);

#	$refseq_id = $db_entry->display_id();

#	print join("\t", $db_entry->display_id() , $db_entry->dbname() ) ."\n";
	
	my $e_id = 1;
	
	foreach my $exon ( @{ $transcript->get_all_translateable_Exons() } ) {
	  

	  $seq_region = $exon->slice->seq_region_name();
	  my $start      = $exon->start();
	  my $end        = $exon->end();
	  ($start, $end) = ($end, $start) if ( $start > $end );
	  $strand     = $exon->strand();
	  
	  $exons{ $refseq_id }{$start} = $end;
	  $transcript_lengths{ $refseq_id } += $end - $start + 1;
#	  print " $refseq_id -- $transcript_lengths{ $refseq_id } $e_id\n";
	  $e_id++;
	}
      }
    }

    next if ( ! $seq_region );


#    print Dumper( \%transcript_lengths );
#    print Dumper( \%exons );
    
    my %used_exons;

    foreach my $refseq_id ( sort { $transcript_lengths{ $b } <=> $transcript_lengths{ $a } } keys %transcript_lengths) {
#      print "$transcript -> $transcript_lengths{ $transcript } \n";

      my $counter = 1;
      $counter = keys %{$exons{$refseq_id}}  if ( $strand == -1 );
      foreach my $start ( sort {$a <=> $b } keys %{$exons{$refseq_id}} ) {
	my $end = $exons{$refseq_id}{$start};

	my $ccds_id = $refseq2ccds{ $refseq_id };

	if ( ! $used_exons{ $start }{ $end }) {
	  print join("\t", $seq_region, $start, $end , $gene_name ."_Exon". $counter, $refseq_id, $ccds_id) . "\n";
	  $used_exons{ $start}{$end}++;
	}

	$counter++ if ( $strand == 1 );
	$counter-- if ( $strand == -1 );
      }
    }

#    print Dumper( \%used_exons );

#    exit;

    # my $counter = 1;
    # $counter = keys %exons  if ( $strand == -1 );
    # foreach my $start ( sort  { $a <=> $b } keys %exons) {
    #   print join("\t", $seq_region, $start, $exons{$start} , $gene_name ."_Exon". $counter, $refseq_id, $ccds_id) . "\n";
    #   $counter++ if ( $strand == 1 );
    #   $counter-- if ( $strand == -1 );
    # }
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



