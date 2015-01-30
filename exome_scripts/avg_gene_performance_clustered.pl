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


my $gene_name = shift;



if ( $gene_name ) {
  print_gene_performance( $gene_name);
}
else {
  my %genes;
  map {$$_{exon_name} =~ s/_Exon.*//ig; $genes{ $$_{ 'exon_name' }}++ } CTRU::depthDB::fetch_exon_names();

#  print Dumper( \%genes );

  foreach my $gene ( sort keys %genes ) {
    print_gene_performance( $gene);
  }

}



# 
# 
# 
# Kim Brugger (09 Dec 2013)
sub print_gene_performance {
  my  ( $gene_name ) = @_;
  
  my @exons = CTRU::depthDB::fetch_regions_by_gene($gene_name);

#  print Dumper( \@exons );
#  exit;


  foreach my $exon ( sort{ my $A = $$a{exon_name};
			   my $B = $$b{exon_name};
			   $A =~ s/.*exon//i;
			   $B =~ s/.*exon//i;
			   $A <=> $B  &&
			   $$a{refseq} cmp $$b{refseq} &&
			   $$a{CCDS} cmp $$b{CCDS} } @exons
      ) {

    my $refseq = $$exon{ refseq };
    my $CCDS   = $$exon{ CCDS   };
    my $rid    = $$exon{ rid   };

    my $transcript_length = 0;
    my $transcript_missed = 0;

    my $region = CTRU::depthDB::fetch_region_hash( $rid );

    my @coverages = CTRU::depthDB::fetch_coverages_by_rid( $rid );

    my (@reads, @min_depth, @missing, @depth_1to5, @depth_6to9, @depth_10to19);
    foreach my $coverage ( @coverages ) {
      my $sample_sequence_hash = CTRU::depthDB::fetch_sample_hash( $$coverage{ sid } );
      next if ( ! $$sample_sequence_hash{ total_reads });
      print join("\t", $$exon{name}, 
		 $$sample_sequence_hash{ name },
		 $$sample_sequence_hash{ total_reads }, 
		 $$sample_sequence_hash{ mapped_reads } - $$sample_sequence_hash{ duplicate_reads }, 
		 $$coverage{min}) . "\n" if ( 0 );

      $transcript_length = $$region{end} - $$region{start}+1;
      push @min_depth, $$coverage{min} || 0;
      push @missing,      range2length($$coverage{missing}) || 0;
      push @depth_1to5,   range2length($$coverage{'1to5'}) || 0;
      push @depth_6to9,   range2length($$coverage{'6to9'}) || 0;
      push @depth_10to19, range2length($$coverage{'10to19'}) || 0;

      $transcript_missed += range2length($$coverage{missing}) + range2length($$coverage{'1to5'}) + range2length($$coverage{'6to9'})+range2length($$coverage{'10to19'});
      push @reads, $$sample_sequence_hash{ total_reads };
#    push @reads, $$sample_sequence_hash{ mapped_reads } - $$sample_sequence_hash{ duplicate_reads };
    }
    
#  print "$$exon{'name'}:::  @min_depth \n";

    
    if (! @min_depth) {
      printf("$$exon{name}\t0\n");
    }
    else {

      my $exp_depth_30M = 0;
      my $reads_for_20x = 0;
      if ( 0 ) {
	my( $intercept, $slope) = linear_regression(\@reads, \@min_depth);
      
	$exp_depth_30M = $intercept + $slope*30_000_000;
	$reads_for_20x = (20 - $intercept) / $slope if ( $slope > 0);
	$reads_for_20x = 0 if ( $slope <= 0);

      
	$reads_for_20x /= 1000000;
	printf("$gene_name\t$refseq\t$CCDS\t$$region{chr}:$$region{start}-$$region{end}\t$$exon{exon_name}\t%.2f\t%dM\n", $exp_depth_30M, $reads_for_20x);
      }
      else {
	my ($meanX, $meanY) = means( \@reads, \@min_depth);
	$exp_depth_30M = sprintf("%.4f", $meanY/$meanX*30_000_000);

	my ($meanXmissing, $meanYmissing) = means( \@reads, \@missing);
	my $exp_missing_30M = sprintf("%.4f", $meanYmissing/$meanXmissing*30_000_000);

	my ($meanX1to5, $meanY1to5) = means( \@reads, \@depth_1to5);
	my $exp_1to5_30M = sprintf("%.4f", $meanY1to5/$meanX1to5*30_000_000);

	my ($meanX6to9, $meanY6to9) = means( \@reads, \@depth_6to9);
	my $exp_6to9_30M = sprintf("%.4f", $meanY6to9/$meanX6to9*30_000_000);

	my ($meanX10to19, $meanY10to19) = means( \@reads, \@depth_10to19);
	my $exp_10to19_30M = sprintf("%.4f", $meanY10to19/$meanX10to19*30_000_000);
  
	my $sd = 0;
	#print "$sample_min_depth\n"; 
	($exp_depth_30M, $sd) = normalise(\@reads, \@min_depth, 30_000_000);
	
	printf("$gene_name\t$refseq\t$CCDS\t$$region{chr}:$$region{start}-$$region{end}\t$$exon{exon_name}\t%.2f\t%dM\t%.2f%%\t$exp_missing_30M\t$exp_1to5_30M\t$exp_6to9_30M\t$exp_10to19_30M\n", $exp_depth_30M, $reads_for_20x, 100*($exp_missing_30M+$exp_1to5_30M+$exp_6to9_30M+$exp_10to19_30M)/$transcript_length);
      }


#      printf("$$exon{chr}:$$exon{start}-$$exon{end}\t$$exon{name}\t%.2f\t%dM\n", $exp_depth_30M, $reads_for_20x);
    }
    
#  print Dumper( \@coverages );
#  exit;
  }
}



# 
# 
# 
# Kim Brugger (24 Mar 2014)
sub range2length {
  my ($ranges ) = @_;

  my $length = 0;

  foreach my $range (split(",", $ranges)) {
    $range =~ s/.*?://;
    my ( $start, $end) = split("-", $range);
    $length += $end- $start + 1;
  }

  return $length;
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

#  print "  my $slope = ($N*$sumXY - ($sumX*$sumY))/( $N*$sumXX - $sumX * $sumX);\n";
  
  my $intercept = ($sumY - $slope*$sumX) / $N;

#  print "$slope $intercept\n";

  return ($intercept, $slope);
  
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
