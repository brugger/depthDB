package CTRU::depthDB;
# 
# 
# 
# 
# Kim Brugger (19 Dec 2013), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;
use POSIX qw( strftime );


use EASIH::DB;

my $dbi;

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub connect {
  my ($dbname, $dbhost, $db_user, $db_pass) = @_;
  $dbhost  ||= "mgsrv01";
  $db_user ||= 'easih_ro';

  $dbi = EASIH::DB::connect($dbname,$dbhost, $db_user, $db_pass);
  return $dbi;
}


#================== Sample functions =========================

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub add_sample {
  my ($name) = @_;

  if ( ! $name ) { 
    print STDERR "add_sample: No sequence name provided\n";
    return -1;
  }

  
  my $sid = fetch_sample_id( $name );
  return $sid if ( $sid );
     
  my %call_hash = ( name => $name);

  return (EASIH::DB::insert($dbi, "sample", \%call_hash));
}



# 
# 
# 
# Kim Brugger (3 Dec 2013)
sub add_sample_stats {
  my ($sid, $total_reads, $mapped_reads, $duplicate_reads, $mean_isize, $gender ) = @_;

  my $ss_name = fetch_sample_name($sid);
  if ( ! $ss_name  ) {
    print STDERR "add_sample_variant: Unknown sample-id $sid '$ss_name'\n";
    return -6;
  }
     
  my %call_hash = ( sid => $sid);
  $call_hash{ total_reads }     = $total_reads     if ( $total_reads     );
  $call_hash{ mapped_reads }    = $mapped_reads    if ( $mapped_reads    );
  $call_hash{ duplicate_reads } = $duplicate_reads if ( $duplicate_reads );
  $call_hash{ mean_isize }      = $mean_isize      if ( $mean_isize      );
  $call_hash{ gender }          = $gender          if ( $gender          );

  return (EASIH::DB::update($dbi, "sample", \%call_hash, "sid"));
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_sample_id {
  my ( $name ) = @_;
  if ( ! $name ) { 
    print STDERR "fetch_sample_id: No sample name provided\n";
    return -1;
  }
  my $q    = "SELECT sid FROM sample WHERE name = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  my @line = EASIH::DB::fetch_array( $dbi, $sth, $name );
  return $line[0] || undef;
}

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_sample_name {
  my ( $sid ) = @_;

  if ( ! $sid ) { 
    print STDERR "fetch_sample_name: No sample id provided\n";
    return "";
  }

  my $q    = "SELECT name FROM sample WHERE sid = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  my @line = EASIH::DB::fetch_array( $dbi, $sth, $sid );
  return $line[0] || undef;
}

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_sample_hash {
  my ( $sid ) = @_;
  if ( ! $sid ) { 
    print STDERR "fetch_sample_hash: No sample id provided\n";
    return {};
  }
  my $q    = "SELECT * FROM sample WHERE sid = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return( EASIH::DB::fetch_hash( $dbi, $sth, $sid ));
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_samples {
  my $q    = "SELECT * FROM sample";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return( EASIH::DB::fetch_array_hash( $dbi, $sth ));
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub update_sample {
  my ($sid, $name, $total_reads, $mapped_reads, $duplicate_reads, $mean_isize, ) = @_;


  if ( ! $sid ) { 
    print STDERR "update_sample: No sample sequence id provided\n";
    return -1;
  }

  if ( ! $name ) { 
    print STDERR "update_sample: No name provided\n";
    return -1;
  }

  my %call_hash;
  $call_hash{sid}               = $sid;
  $call_hash{name}              = $name            if ($name);
  $call_hash{ total_reads }     = $total_reads     if ( $total_reads );
  $call_hash{ mapped_reads }    = $mapped_reads    if ( $mapped_reads );
  $call_hash{ duplicate_reads } = $duplicate_reads if ( $duplicate_reads );
  $call_hash{ mean_isize }      = $mean_isize      if ( $mean_isize );
  return (EASIH::DB::update($dbi, "sample", \%call_hash, "sid"));
}


#================== region functions =========================

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub add_region {
  my ($chr, $start, $end) = @_;

  if ( ! $chr ) { 
    print STDERR "add_region: No chr provided\n";
    return -1;
  }

  if ( ! $start ) { 
    print STDERR "add_region: No region start provided\n";
    return -2;
  }

  if ( ! $end ) { 
    print STDERR "add_region: No region end provided\n";
    return -3;
  }

  my $rid = fetch_region_id_by_position( $chr, $start, $end );
  return $rid if ( $rid );
     
  my %call_hash = ( chr       => $chr,
		    start     => $start,
		    end       => $end);

  return (EASIH::DB::insert($dbi, "region", \%call_hash));
}

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_region_hash {
  my ( $rid ) = @_;
  if ( ! $rid ) { 
    print STDERR "fetch_region_hash: No region id provided\n";
    return {};
  }
  my $q    = "SELECT * FROM region WHERE rid = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return( EASIH::DB::fetch_hash( $dbi, $sth, $rid ));
}



# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_regions_by_gene {
  my ( $gene ) = @_;

  if ( ! $gene ) { 
    print STDERR "fetch_region_by_gene: No gene name provided\n";
    return -1;
  }
  my $q    = "SELECT * FROM transcript WHERE exon_name like ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return EASIH::DB::fetch_array_hash( $dbi, $sth, "$gene\_E%" );
}




# 
# 
# 
# Kim Brugger (19 Dec 2013)
sub fetch_exon_names {

  my $q    = "SELECT exon_name FROM transcript";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return EASIH::DB::fetch_array_hash( $dbi, $sth);
  
}

# 
# 
# 
# Kim Brugger (19 Dec 2013)
sub fetch_regions {

  my $q    = "SELECT * FROM region";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return EASIH::DB::fetch_array_hash( $dbi, $sth);
  
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_region_id_by_position {
  my ( $chr, $start, $end ) = @_;

  if ( ! $chr || !$start || !$end) { 
    print STDERR "fetch_region_id_by_position: No chromosome or position provided\n";
    return -1;
  }
  my $q    = "SELECT rid FROM region WHERE chr = ? AND start = ? and end  = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  my @line = EASIH::DB::fetch_array( $dbi, $sth, $chr, $start, $end );
  return $line[0] || undef;
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub update_region {
  my ($rid, $chr, $start, $end) = @_;

  if ( ! $rid ) { 
    print STDERR "update_region: No sample sequence id provided\n";
    return -1;
  }

  my %call_hash;
  $call_hash{rid}       = $rid   if ($rid);
  $call_hash{chr}       = $chr   if ($chr);
  $call_hash{start}     = $start if ($start);
  $call_hash{end}       = $end   if ($end);

  return (EASIH::DB::update($dbi, "region", \%call_hash, "rid"));
}


#================== coverage functions =========================

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub add_coverage {
  my ($sid, $rid, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19) = @_;

  if ( ! $sid ) { 
    print STDERR "add_coverage: No sample id provided\n";
    return -1;
  }

  if ( ! $rid ) { 
    print STDERR "add_coverage: No region id provided\n";
    return -2;
  }

  if ( ! defined $min ) { 
    print STDERR "add_coverage: No min depth provided\n";
    return -3;
  }

  if ( ! defined  $mean ) { 
    print STDERR "add_coverage: No mean depth provided\n";
    return -4;
  }

  if ( ! defined  $max ) { 
    print STDERR "add_coverage: No max depth provided\n";
    return -5;
  }

  my $ss_name = fetch_sample_name( $sid );
  if ( ! $ss_name  ) {
    print STDERR "add_coverage: Unknown sample-id: $sid\n";
    return -8;
  }

  my $r_hash = fetch_region_hash( $rid );
  if ( ! $r_hash || keys %{$r_hash} == 0) {
    print STDERR "add_coverage: Unknown region-id: $rid\n";
    return -9;
  }
  
  my $c_hash = fetch_coverage_hash( $sid, $rid );
  return 1 if ( $c_hash && keys %{$c_hash} > 0 );

  my %call_hash = ( sid     => $sid,
		    rid     => $rid,
		    min     => $min,
		    mean    => $mean,
		    max     => $max);


  return (EASIH::DB::insert($dbi, "coverage", \%call_hash));
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub add_coverages {
  my ( @coverages ) = @_;

  if ( ! @coverages ) { 
    print STDERR "add_coverages: No coverage(s) provided\n";
    return -1;
  }

  my @new_coverages;
  foreach my $entry ( @coverages ) {
    my ($sid, $rid, $min, $mean, $max, $missing, $depth_1to5, $depth_6to9, $depth_10to19) = @$entry;



    if ( ! $sid ) { 
      print STDERR "add_coverage: No sample id provided\n";
      return -1;
    }

    if ( ! $rid ) { 
      print STDERR "add_coverage: No region id provided\n";
      return -2;
    }


    my $ss_name = fetch_sample_name( $sid );
    if ( ! $ss_name  ) {
      print STDERR "add_coverage: Unknown sample-id: $sid\n";
      return -8;
    }

    my $r_hash = fetch_region_hash( $rid );
    if ( ! $r_hash || keys %{$r_hash} == 0) {
      print STDERR "add_coverage: Unknown region-id: $rid\n";
      return -9;
    }
  
    my $c_hash = fetch_coverage_hash( $sid, $rid );
    return 1 if ( $c_hash && keys %{$c_hash} > 0 );
    
    my %call_hash = ( 'sid'     => $sid,
		      'rid'     => $rid,
		      'min'     => $min,
		      'mean'    => $mean,
		      'max'     => $max,
		      'missing' => $missing,
		      '1to5'    => $depth_1to5, 
		      '6to9'    => $depth_6to9,  
		      '10to19'  => $depth_10to19);
    push @new_coverages, \%call_hash;

    if ( int( @new_coverages ) == 500 ) {
      (EASIH::DB::insert($dbi, "coverage", \@new_coverages));
      undef @new_coverages;
    }
  }

  if ( int( @new_coverages ) ) {
    (EASIH::DB::insert($dbi, "coverage", \@new_coverages));
  }

  return;
}



# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_coverage_hash {
  my ($sid,  $rid ) = @_;
  if ( ! $rid || ! $sid ) { 
    print STDERR "fetch_coverage_hash: No variant and/or sample id prorided\n";
    return {};
  }
  my $q    = "SELECT * FROM coverage WHERE sid = ? AND rid = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return( EASIH::DB::fetch_hash( $dbi, $sth, $sid, $rid ));
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_coverages_by_rid {
  my ($rid ) = @_;
  if ( ! $rid ) { 
    print STDERR "fetch_coverages_by_rid: No region id prorided\n";
    return {};
  }

  my $q    = "SELECT * FROM coverage WHERE rid = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return( EASIH::DB::fetch_array_hash( $dbi, $sth, $rid ));
}


# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub update_coverage {
  my ($sid, $rid, $min, $mean, $max, $lows, $missing) = @_;

  if ( ! $sid ) { 
    print STDERR "add_coverage: No sequence id provided\n";
    return -1;
  }

  if ( ! $rid ) { 
    print STDERR "add_coverage: No region id provided\n";
    return -2;
  }

  my $ss_name = fetch_sample_name( $sid );
  if ( ! $ss_name  ) {
    print STDERR "add_coverage: Unknown sample-id: $sid\n";
    return -3;
  }

  my $r_hash = fetch_region_hash( $rid );
  if ( ! $r_hash || keys %{$r_hash} == 0) {
    print STDERR "add_coverage: Unknown region-id: $rid\n";
    return -4;
  }
  
  my %call_hash;
  $call_hash{ sid }      = $sid     if ( $sid             );
  $call_hash{ rid  }     = $rid     if ( $rid             );
  $call_hash{ min  }     = $min     if ( $min             );
  $call_hash{ mean }     = $mean    if ( $mean            );
  $call_hash{ max }      = $max     if ( $max             );
  $call_hash{ lows  }    = $lows    if ( defined $lows    );
  $call_hash{ missing }  = $missing if ( defined $missing );

  return (EASIH::DB::update($dbi, "coverage", \%call_hash, "sid", "rid"));
}


#================== transcript functions =========================

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub add_transcript {
  my ($rid, $exon_name, $refseq, $CCDS) = @_;

  if ( ! $rid ) { 
    print STDERR "add_transcript: No region id provided\n";
    return -1;
  }

  my $r_hash = fetch_region_hash( $rid );
  if ( ! $r_hash || keys %{$r_hash} == 0) {
    print STDERR "add_transcript: Unknown region-id: $rid\n";
    return -9;
  }
  
#  my $r_hash = fetch_region_hash( $rid );
#  return 1 if ( $c_hash && keys %{$c_hash} > 0 );

  my $transcript_hash = fetch_transcript_hash( $rid, $exon_name, $refseq, $CCDS );
  return 1 if ( $transcript_hash );
#  die Dumper( $transcript_hash );

  my %call_hash = ( rid       => $rid,
		    exon_name => $exon_name,
		    refseq    => $refseq,
		    CCDS      => $CCDS);

  return (EASIH::DB::insert($dbi, "transcript", \%call_hash));
}

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub fetch_transcript_hash {
  my ($rid, $exon_name, $refseq, $CCDS ) = @_;
  if ( ! $rid ) { 
    print STDERR "fetch_transcript_hash: No region id prorided\n";
    return {};
  }
  my $q    = "SELECT * FROM transcript WHERE rid = ? and exon_name = ? and refseq = ? and CCDS = ?";
  my $sth  = EASIH::DB::prepare($dbi, $q);
  return( EASIH::DB::fetch_hash( $dbi, $sth, $rid, $exon_name, $refseq, $CCDS ));
}

# 
# 
# 
# Kim Brugger (20 Nov 2013)
sub update_transcript {
  my ($rid, $exon_name, $refseq, $CCDS) = @_;

  if ( ! $rid ) { 
    print STDERR "add_transcript: No region id provided\n";
    return -1;
  }

  my $r_hash = fetch_region_hash( $rid );
  if ( ! $r_hash || keys %{$r_hash} == 0) {
    print STDERR "add_transcript: Unknown region-id: $rid\n";
    return -9;
  }
  
#  my $r_hash = fetch_region_hash( $rid );
#  return 1 if ( $c_hash && keys %{$c_hash} > 0 );


  my %call_hash;
  $call_hash{ rid  }      = $rid       if ( $rid             );
  $call_hash{ exon_name } = $exon_name if ( $exon_name             );
  $call_hash{ refseq }     = $refseq   if ( $refseq            );
  $call_hash{ CCDS }      = $CCDS      if ( $CCDS             );

  return (EASIH::DB::update($dbi, "transcript", \%call_hash, "rid"));
}




1;
