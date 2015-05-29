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

use Spreadsheet::WriteExcel;


my $MIN_READS = 25_000_000;

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
#$dbname = "frags_exome";

my $dbi = CTRU::depthDB::connect($dbname, $dbhost, "easih_admin", "easih");

my $gene_name = shift;
my $sample_name = $opts{s};

my $excel_file = "$sample_name.dosage.xls";

my $workbook = Spreadsheet::WriteExcel->new( $excel_file );

my $bold     = $workbook->add_format(bold => 1);
my $red_cell    = $workbook->add_format(color => 'red'    );
my $orange_cell = $workbook->add_format(color => 'orange' );
my $green_cell  = $workbook->add_format(color => 'green'  );
my $format_top_left    = $workbook->add_format(top => 1, left=>1);
my $format_left        = $workbook->add_format(left=>1);
my $format_left_bottom = $workbook->add_format(left=>1, bottom=>1);
my $format_bottom = $workbook->add_format(bottom=>1);

my $dosage_tab   = $workbook->add_worksheet('dosage');

$dosage_tab->write(0, 0, "Sample", $format_bottom);
$dosage_tab->write(0, 1, "$sample_name", $format_bottom);

my $dosage_offset = 2;

my @fields = ("Exon name", "Effect", "Confidence");
for(my$i=0;$i< @fields; $i++ ) {
  $dosage_tab->write($dosage_offset, $i, $fields[$i], $bold);
}
$dosage_offset++;

foreach my $gene_name ( keys %gene_list ) {
  
  my @exons = CTRU::depthDB::fetch_regions_by_gene($gene_name);



  
  foreach my $exon ( sort{ my $A = $$a{exon_name};
			   my $B = $$b{exon_name};
			   $A =~ s/.*exon//i;
			   $B =~ s/.*exon//i;
			   $A <=> $B  
		     } @exons
      ) {
  
    my $refseq = $$exon{ refseq };
    my $CCDS   = $$exon{ CCDS   };
    my $rid    = $$exon{ rid   };

    next if ( $refseq !~ /$gene_list{ $gene_name }/);

    my $region = CTRU::depthDB::fetch_region_hash( $rid );

    my $region_size = $$region{ end } - $$region{ start } + 1;
  
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
    

    next if (! @min_depth);
    if (!$sample_reads) {
      print "No sample reads for $sample_name\n";
      exit;
    }
    my ($exp_depth, $sd) = normalise(\@reads, \@min_depth, $sample_reads);
    next if ( $exp_depth == 0);


    my $ratio = sprintf("%.4f", $sample_min_depth/$exp_depth);
    $ratio = 0.0001 if ( $ratio == 0);
    my $log_ratio = sprintf("%.4f", log( $ratio )/log(2));
    
#    printf("$$exon{ exon_name }\t$sample_min_depth\t%.2f\t%.2f\t$log_ratio\t$ratio\n", $exp_depth, $sd);
    my $mean = $exp_depth;

    if ( $sample_min_depth < $mean-2.58*$sd ) {


      if ( $ratio < 0.15 ) {
	print "$$exon{ exon_name } homo deletion ( > 99% confidence) RATIO:$ratio $region_size\n";
	$dosage_tab->write($dosage_offset, 1, "deletion (homo)");
      }
      else {
	print "$$exon{ exon_name } het deletion ( > 99% confidence) RATIO:$ratio $region_size\n";
	$dosage_tab->write($dosage_offset, 1, "deletion (het)");
      }

      $dosage_tab->write($dosage_offset, 0, "$$exon{ exon_name }");
      $dosage_tab->write($dosage_offset, 2, "> 99%");
      $dosage_offset++;


    }
    elsif ( $sample_min_depth < $mean-1.96*$sd ) {
      if ( $ratio < 0.15 ) {
	print "$$exon{ exon_name } homo deletion ( > 95% confidence) RATIO:$ratio $region_size\n";
	$dosage_tab->write($dosage_offset, 1, "deletion (homo)");
      }
      else {
	print "$$exon{ exon_name } het deletion ( > 95% confidence) RATIO:$ratio $region_size\n";
	$dosage_tab->write($dosage_offset, 1, "deletion (het)");
      }
      $dosage_tab->write($dosage_offset, 0, "$$exon{ exon_name }");
      $dosage_tab->write($dosage_offset, 2, "> 95%");
      $dosage_offset++;
    }
    elsif  ( $sample_min_depth > $mean+2.58*$sd ) {
      print "$$exon{ exon_name } duplication ( > 99% confidence) RATIO:$ratio $region_size\n";
      $dosage_tab->write($dosage_offset, 0, "$$exon{ exon_name }");
      $dosage_tab->write($dosage_offset, 1, "Duplication");
      $dosage_tab->write($dosage_offset, 2, "> 95%");
      $dosage_offset++;
    }
    elsif ( $sample_min_depth > $mean+1.96*$sd ) {
      print "$$exon{ exon_name } duplication ( > 95% confidence) RATIO:$ratio $region_size\n";
      $dosage_tab->write($dosage_offset, 0, "$$exon{ exon_name }");
      $dosage_tab->write($dosage_offset, 1, "duplication");
      $dosage_tab->write($dosage_offset, 2, "> 99%");
      $dosage_offset++;
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

