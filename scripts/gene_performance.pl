#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (11 May 2015), contact: kbr@brugger.dk

use strict;
use warnings;
use Data::Dumper;

use lib '/software/packages/depthDB/modules';
use CTRU::depthDB;

my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';
my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");

my $MIN_READS  = 25_000_000;
my $NORM_READS = 30_000_000;


while( my $gene_name = <> ) {
  if ($gene_name =~ / / || $gene_name =~ /^\Z/) {
    print $gene_name;
    next;
  }
  chomp $gene_name;
  my %expected_coverages;

  next if ($gene_name =~ /^#/);

  my @exons = CTRU::depthDB::fetch_regions_by_gene($gene_name);

  foreach my $exon ( @exons ) {
      
    my $refseq = $$exon{ refseq };
    my $CCDS   = $$exon{ CCDS   };
    my $rid    = $$exon{ rid   };

  
    my $region = CTRU::depthDB::fetch_region_hash( $rid );
    
    my @coverages = CTRU::depthDB::fetch_coverages_by_rid( $rid );
    
    my (@reads, @min_depth);
    my (@plotting);
    my @exon_performance;
    
    foreach my $coverage ( @coverages ) {
      my $sample_sequence_hash = CTRU::depthDB::fetch_sample_hash( $$coverage{ sid } );
      next if ( $$sample_sequence_hash{ name } =~ /beta/i);

      next if ( ! $$sample_sequence_hash{ total_reads }|| $$sample_sequence_hash{ total_reads } < $MIN_READS);


      my $low_or_missing = range2length($$coverage{'missing'});
      $low_or_missing +=  range2length($$coverage{'1to5'});  
      $low_or_missing +=  range2length($$coverage{'6to9'});  
      $low_or_missing +=  range2length($$coverage{'10to19'});
      my $region_length = $$region{ 'end' }-$$region{ 'start' } + 1;

      my $well_covered = ($region_length-$low_or_missing)*100/$region_length;
#      print "$well_covered $low_or_missing*100/$region_length\n";
#      exit if ( ! $well_covered );
      push @exon_performance, $well_covered;
#      push @min_depth, $$coverage{mean};
      $$sample_sequence_hash{ duplicate_reads } ||= 0;
      push @reads, $$sample_sequence_hash{ mapped_reads } - $$sample_sequence_hash{ duplicate_reads };
    }  
    

    push @{$expected_coverages{ $$exon{ 'refseq'}}}, mean_exp_coverage( @exon_performance );
  }
  foreach my $transcript ( sort keys %expected_coverages ) {
    print join("\t", $gene_name, $transcript, mean_exp_coverage( @{$expected_coverages{ $transcript}}))."\n";
  }
}





# 
# 
# 
# Kim Brugger (11 May 2015)
sub mean_exp_coverage {
  my (@values ) = @_;
  my $sum = 0;
  map { $sum += $_ } @values;

  return sprintf("%.2f", $sum/int(@values));
  
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
# Kim Brugger (24 Mar 2014)
sub range2length {
  my ($ranges ) = @_;

  return 0 if (! $ranges || $ranges eq "");

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

#  print "  my $slope = ($N*$sumXY - ($sumX*$sumY))/( $N*$sumXX - $sumX*$sumX);\n";
  
  my $intercept = ($sumY - $slope*$sumX) / $N;

#  print "$intercept = ($sumY - $slope*$sumX) / $N;\n";

#  print "$slope $intercept\n";

  return ($intercept, $slope);
  
}

