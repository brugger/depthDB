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

my $genes2transcripts_file = "/data/gemini/genes2transcripts";

my %genes2transcripts;
open(my $in, $genes2transcripts_file);
map{ my @F = split(/\s+/, $_); $genes2transcripts{ uc($F[0])} = uc($F[1]) if ( $F[1] )} <$in>;
close($in);

my $genepanel_file         = "/data/gemini/genepanels";
my %genepanels;
open( $in, $genepanel_file) || die "Could not open '$genepanel_file': $!\n";
while( <$in> ) {
  chomp;
  my @F = split(/\t/, $_); 
  push @{$genepanels{ uc($F[0])}}, uc($F[3]) if ( $F[0] && $F[3]);
  push @{$genepanels{ uc($F[1])}}, uc($F[3]) if ( $F[1] && $F[3]);
  push @{$genepanels{ uc($F[2])}}, uc($F[3]) if ( $F[2] && $F[3]);
}
close($in);


use Getopt::Std;
my $opts = 's:G:';
my %opts;
getopts($opts, \%opts);

my %gene_list;
%gene_list = readin_gene_list($opts{ 'G' });

#exit;
use Getopt::Std;
my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';

my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");
foreach my $gene_name ( keys %gene_list ) {
  
  my @exons = CTRU::depthDB::fetch_regions_by_gene($gene_name);



  
  foreach my $exon ( sort{ my $A = $$a{exon_name};
			   my $B = $$b{exon_name};
			   $A =~ s/.*exon//i;
			   $B =~ s/.*exon//i;
			   $A <=> $B  
		     } @exons
      ) {

    my $rid    = $$exon{ rid   };
    my $region = CTRU::depthDB::fetch_region_hash( $rid );

    print join("\t", $$region{chr}, $$region{ start}, $$region{ end }, $$region{ rid }). "\n";
  }
}


# 
# 
# 
# Kim Brugger (27 Aug 2013)
sub readin_gene_list {
  my ($gene_list) = @_;

  my %gene_list = ();

  open( my $in, $gene_list ) || die "Could not open '$gene_list': $!\n";
  while(<$in>) {
    chomp;
    my ( $gene, $transcript) = split("\t");
    next if ( /^\z/ );

#    print Dumper( \%genepanels );

    if ( $genepanels{ uc( $gene ) } ) {
      foreach my $panel_gene ( @{$genepanels{ uc( $gene ) }} ) {
	my $panel_transcript ||= $genes2transcripts{uc($panel_gene)};
	die "No panel-transcript for '$panel_gene' \n" if ( ! $panel_transcript || $panel_transcript eq "");
	$panel_transcript ||= "-";
	
	$gene_list{ uc($panel_gene) } = uc($panel_transcript);
	
      }
      next;
    }

    $transcript ||= $genes2transcripts{uc($gene)};

    die "No transcript for $gene \n" if ( ! $transcript || $transcript eq "");
    
    $gene_list{ uc($gene) } = uc($transcript);
  }
  return %gene_list;
}


# 
# 
# 
# Kim Brugger (21 Feb 2014)
sub means {
  my ( $xs, $ys ) = @_;

  my $sumX = 0;
  my $sumY = 0;

  map { $sumX += $_ } @$xs;
  map { $sumY += $_ } @$ys;

  return ( $sumX/@$xs, $sumY/@$ys)
  
}

sub print_depths {
   my ( $xs, $ys) = @_;

  my @normals;
  my $summed_normals = 0;

  for( my $i = 0; $i < @$xs;$i++) {
    print join("\t", $$xs[ $i ], $$ys[ $i ])."\n";
  }

}



sub normalise {
   my ( $xs, $ys, $normal_reads) = @_;

  my @normals;
  my $summed_normals = 0;

  for( my $i = 0; $i < @$xs;$i++) {
    my $reads = $$xs[ $i ];
    my $depth = $$ys[ $i ];

    my $normal = $depth*$normal_reads/$reads;
    push @normals, $normal;
  }

#  print join("\n", @normals) . "\n";


  for( my $i = 0; $i < @normals; $i++) {
    $summed_normals += $normals[ $i ];
  } 

  return sd( \@normals);

  return sprintf("%.2f", $summed_normals/(int(@$xs)));

#  print join("\n", @normals ) . "\n";

}

sub sd {
  my ( $values ) = @_;

  my $sum  = 0;
  my $sum2 = 0;
  my $value_count = 0;

  foreach my $value ( @$values ) { 

    $sum  += $value; 
    $sum2 += $value*$value;
    $value_count++;
  }
    my $mean = $sum/$value_count;
    my $sd   = sqrt(( $value_count * $sum2-$sum*$sum)/($value_count *( $value_count -1)));

#    print "$sum\t$sum2\t$input_files\n";

   return ($mean, $sd);
}





# 
# 
# 
# Kim Brugger (06 Dec 2013)
sub linear_regression {
  my ( $xs, $ys ) = @_;

  my $N = @$xs;

  my ( $sumX, $sumY, $sumXY, $sumXX) = (0,0,0,0);

  for(my $i=0;$i < $N; $i++ ) {

    $sumX  += $$xs[ $i ];
    $sumY  += $$ys[ $i ];
    $sumXX += $$xs[ $i ] * $$xs[ $i ];
    $sumXY += $$xs[ $i ] * $$ys[ $i ];

  }

#  print join("\t", $sumX, $sumY, $sumXX, $sumXY, $N) . "\n";

  my $slope = ($N*$sumXY - ($sumX*$sumY))/( $N*$sumXX - $sumX*$sumX);

#  print "  my $slope = ($N*$sumXY - ($sumX*$sumY))/( $N*$sumXX - $sumX*$sumX);\n";
  
  my $intercept = ($sumY - $slope*$sumX) / $N;

#  print "$intercept = ($sumY - $slope*$sumX) / $N;\n";

#  print "$slope $intercept\n";

  return ($intercept, $slope);
  
}

