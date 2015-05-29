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


my $MIN_READS = 25_000_000;


#exit;
use Getopt::Std;
my $dbhost = 'mgsrv01';
my $dbname = 'depths_exome_5bp';
#$dbname = "frags_exome";

my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");


my $gene_name = shift;
my $sample_name = shift;


  
my @exons = CTRU::depthDB::fetch_regions_by_gene($gene_name);

#  print Dumper( \@exons );
#  exit;


printf("Exon_name\tSample mean depth\tExp depth\tRatio\tLog Ratio\n");
#printf("Exon_name\tLog_Ratio\n");



foreach my $exon ( sort{ my $A = $$a{exon_name};
			 my $B = $$b{exon_name};
			 $A =~ s/.*exon//i;
			 $B =~ s/.*exon//i;
			 $A <=> $B  
#   		         $$a{refseq} cmp $$b{refseq} &&
#			 $$a{CCDS} cmp $$b{CCDS} 
		   } @exons
    ) {
  
  my $refseq = $$exon{ refseq };
  my $CCDS   = $$exon{ CCDS   };
  my $rid    = $$exon{ rid   };

#  next if ( $$exon{ exon_name} ne "FBN1_Exon35");
#  next if ( $$exon{ exon_name} ne "FBN1_Exon63");
#  next if ( $$exon{ exon_name} ne "F8_Exon14");
  
  my $region = CTRU::depthDB::fetch_region_hash( $rid );
  
  my @coverages = CTRU::depthDB::fetch_coverages_by_rid( $rid );
  
  my (@reads, @min_depth);
  my (@plotting);
  my ($sample_reads, $sample_min_depth);

  foreach my $coverage ( @coverages ) {
    my $sample_sequence_hash = CTRU::depthDB::fetch_sample_hash( $$coverage{ sid } );

#    print "$$sample_sequence_hash{ name }\n";

    next if ( $$sample_sequence_hash{ name } =~ /beta/i);
    if ( $$sample_sequence_hash{ name } eq $sample_name ) {
      $sample_min_depth = $$coverage{mean} || $$coverage{min} || 0;
#      $sample_reads = $$sample_sequence_hash{ total_reads };
      $sample_reads = $$sample_sequence_hash{ mapped_reads } - $$sample_sequence_hash{ duplicate_reads };
      next;
    }

#    next if ( ! $sample_reads );
      
    next if ( ! $$sample_sequence_hash{ total_reads }|| $$sample_sequence_hash{ total_reads } < $MIN_READS);
    push @min_depth, $$coverage{mean} || $$coverage{min} || 0;
#    push @reads, $$sample_sequence_hash{ total_reads };
#    push @reads, $$sample_sequence_hash{ mapped_reads } - $$sample_sequence_hash{ duplicate_reads };
    $$sample_sequence_hash{ duplicate_reads } ||= 0;
    push @reads, $$sample_sequence_hash{ mapped_reads } - $$sample_sequence_hash{ duplicate_reads };
    push @plotting, "$$sample_sequence_hash{ name }, $$sample_sequence_hash{ total_reads }\t$$coverage{min}";
  }
  
#  print "$$exon{'name'}:::  @min_depth \n";


#  print Dumper( $region );

  if (! @min_depth) {
    printf("$$exon{name}\t0\n");
  }
  else {

    my ($exp_depth, $sd) = normalise(\@reads, \@min_depth, $sample_reads);

    next if ( $exp_depth == 0);


    my $ratio = sprintf("%.4f", $sample_min_depth/$exp_depth);
    $ratio = 0.0001 if ( $ratio == 0);
    my $log_ratio = sprintf("%.4f", log( $ratio )/log(2));
    
    printf("$$exon{ exon_name }\t$sample_min_depth\t%.2f\t%.2f sd\tlog2: $log_ratio\t$ratio\n", , $sd);
    my $mean = $exp_depth;

    if ( $sample_min_depth < $mean-2.58*$sd ) {
      print "$$exon{ exon_name } deletion ( > 99% confidence)\n";
    }
    elsif ( $sample_min_depth < $mean-1.96*$sd ) {
      print "$$exon{ exon_name } deletion ( > 95% confidence)\n";
    }
    elsif  ( $sample_min_depth > $mean+2.58*$sd ) {
      print "$$exon{ exon_name } duplication ( > 99% confidence)\n";
    }
    elsif ( $sample_min_depth > $mean+1.96*$sd ) {
      print "$$exon{ exon_name } duplication ( > 95% confidence)\n";
    }
#  printf("$$exon{ exon_name }\t$sample_min_depth\t%.2f\t%2.f\t%.2f\t$ratio\n", $mean-1.96*$sd, $mean, $mean+1.96*$sd);
#  printf("$$exon{ exon_name }\t$sample_min_depth\t%.2f\t%2.f\t%.2f\t$ratio\n", $mean-2.58*$sd, $mean, $mean+2.58*$sd);
  
#  printf("$$exon{ exon_name }\t$log_ratio\n");
  
    
  }

#  print join("\n", @plotting) ."\n";
#  print "\n$sample_reads\t$sample_mean_depth\n"
  
#  print Dumper( \@coverages );
#  exit;
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

