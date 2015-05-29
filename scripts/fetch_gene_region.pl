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

#exit;
use Getopt::Std;
my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';
my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");


my $FLANK = 500;

while( my $gene_name = <> ) {
  
  chomp( $gene_name );


  my @exons = CTRU::depthDB::fetch_regions_by_gene($gene_name);

  my ($chr, $start, $end) = (undef, 0,0);


  
  foreach my $exon ( sort{ my $A = $$a{exon_name};
			   my $B = $$b{exon_name};
			   $A =~ s/.*exon//i;
			   $B =~ s/.*exon//i;
			   $A <=> $B  
		     } @exons
      ) {

    my $rid    = $$exon{ rid   };
    my $region = CTRU::depthDB::fetch_region_hash( $rid );

    if ( !$start ) {
      ($chr, $start, $end) = ($$region{ chr }, $$region{ start }, $$region{ end });
      next;
    }

    $start = $$region{ start} if ( $start > $$region{ start });
    $end   = $$region{ end  } if ( $end   < $$region{ end  });


  }

  $start -= $FLANK;
  $end   += $FLANK;

  print join("\t", $chr, $start, $end, $gene_name). "\n";
}
