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

use lib '/software/packages/bioperl/';


my $species     = "human";
my $host        = 'mgsrv01.medschl.cam.ac.uk';
my $user        = "easih_ro";

#use lib '/software/packages/ensembl/modules/';
use lib '/software/packages/ensembl_74/modules/';

my $reg = 'Bio::EnsEMBL::Registry';
#$reg->load_registry_from_db(-host => $host,-user => $user);
$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');


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


my $exit_count = 9;
# get adaptors
my $sa  = $reg->get_adaptor($species, 'core', 'slice');
my $ga  = $reg->get_adaptor($species, 'core', 'gene');
my $ta  = $reg->get_adaptor($species, 'core', 'transcript');


while (my $gene_id = <>) {
  
  chomp ($gene_id );
#  print "$transcript_id";

  my @genes = @{ $ga->fetch_all_by_external_name( $gene_id )};
#  my @genes = @{ $ga->fetch_all_by_display_label( $gene_id )};

  if ( ! @genes ) {

    print "'$gene_id' not found in the Ensembl database\n";
  }

  my $gene_name = "";
  foreach my $gene ( @genes ) {
    
    my $seq_region = undef;
    my $strand    = 1;
    my $transcripts = $gene->get_all_Transcripts();
    my %exons;

    my $total_refseq_entries = 0;
    my $total_ccds_entries   = 0;
    my $usable_transcripts   = 0;
    
    
    while ( my $transcript = shift @{$transcripts} ) {
#      print "------ NEW TRANSCRIPT ------\n";

      my @db_entries = (@{$transcript->get_all_DBEntries('RefSeq_dna')}, @{ $transcript->get_all_DBEntries('RefSeq_mRNA')},  @{ $transcript->get_all_DBEntries('CCDS')});
#      @db_entries = (@{$transcript->get_all_DBEntries()});

#      print Dumper( \@db_entries );
      
      $gene_name = $gene->external_name();
#      next if ( $gene_name ne $gene_id);

      my ($refseq_entry, $ccds_entry) = (0,0);
      my $ccds_id = "";
      my @refseq_ids;
      
      foreach my $db_ent ( @db_entries ) {
#	print join("\t", $db_ent->display_id() , $db_ent->dbname() ) ."\n";
	if ( $db_ent->dbname() eq "RefSeq_dna" || $db_ent->dbname() eq "RefSeq_mRNA" ) {
	  $refseq_entry++;
	  $total_refseq_entries++;
	}

	push @refseq_ids, $db_ent->display_id() if ( $db_ent->display_id() && $db_ent->dbname() eq "RefSeq_dna");
	push @refseq_ids, $db_ent->display_id() if ( $db_ent->display_id() && $db_ent->dbname() eq "RefSeq_mRNA");

	if ( $db_ent->dbname() eq "CCDS") {
	  $ccds_entry++;
	  $total_ccds_entries++;
	}
	$ccds_id = $db_ent->display_id() if ( $db_ent->dbname() eq "CCDS");
#	print Dumper( \@refseq_ids );
      }

      my $refseq_id = join(",", @refseq_ids);

      $refseq_id ||= $ccds_id;


      next if ( ! $refseq_id );
#      print "REFSEQ/CCDS = $refseq_id/$ccds_id  $refseq_entry/$ccds_entry\n";

      foreach my $exon ( @{ $transcript->get_all_translateable_Exons() } ) {
	  
	
	$seq_region = $exon->slice->seq_region_name();
	next if ( $seq_region =~ /patch/i);
	my $start      = $exon->start();
	my $end        = $exon->end();
	($start, $end) = ($end, $start) if ( $start > $end );
	$strand     = $exon->strand();
	  
	$exons{ $refseq_id }{$start} = $end;
      }

#      print Dumper( $exons{ $refseq_id } );
      $usable_transcripts++;
      next if ( ! $refseq_entry && ! $ccds_entry);

      my $counter = 1;
      $counter = keys %{$exons{$refseq_id}}  if ( $strand == -1 );
      foreach my $start ( sort {$a <=> $b } keys %{$exons{$refseq_id}} ) {
	my $end = $exons{$refseq_id}{$start};
	
#	my $ccds_id = $refseq2ccds{ $refseq_id };
	
	print join("\t", $seq_region, $start, $end , $gene_name ."_Exon". $counter, $refseq_id, $ccds_id) . "\n";
	  
	$counter++ if ( $strand == 1 );
	$counter-- if ( $strand == -1 );
      }
    }


#    print "! $usable_transcripts && ($total_refseq_entries || $total_ccds_entries\n";

#      print Dumper( \%exons );
    if ( ! $usable_transcripts && ($total_refseq_entries || $total_ccds_entries )) {
#      print Dumper( \%exons );      
      
      foreach my $id ( keys %exons ) {
	my ( $refseq_id, $ccds_id) = ("",'');
	$refseq_id = $id if ( $id =~ /NM/);
	$ccds_id   = $id if ( $id =~ /CCDS/);
	my $counter = 1;
	$counter = keys %{$exons{$id}}  if ( $strand == -1 );
	foreach my $start ( sort {$a <=> $b } keys %{$exons{$id}} ) {
	  my $end = $exons{$id}{$start};
	  
#	my $ccds_id = $refseq2ccds{ $refseq_id };


	  
	  print join("\t", $seq_region, $start, $end , $gene_name ."_Exon". $counter, $refseq_id, $ccds_id) . "\n";
	  
	  $counter++ if ( $strand == 1 );
	  $counter-- if ( $strand == -1 );
	}
      }




    }
  }
}


# 
# 
# 
# Kim Brugger (04 Feb 2014)
sub identical_gene_set {
  my ( $ref1, $ref2 ) = @_;

  foreach my $pos ( keys %$ref1 ) {

    print "$$ref1{ $pos } ne $$ref2{ $pos } \n"   if ( $$ref1{ $pos } ne $$ref2{ $pos } );

    return 0   if ( $$ref1{ $pos } ne $$ref2{ $pos } );
  }

  foreach my $pos ( keys %$ref2 ) {
    print "$$ref1{ $pos } ne $$ref2{ $pos } \n"   if ( $$ref1{ $pos } ne $$ref2{ $pos } );
    return 0   if ( $$ref1{ $pos } ne $$ref2{ $pos } );
  }

  return 1;
  
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



