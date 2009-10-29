package best_blast_hit;

use strict;
use genomedbrange;
use Bio::Root::Root;
use Bio::Tools::BPlite;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'best_blast_hit';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my $coverage_threshold = $self->_rearrange([qw(COVERAGE_THRESHOLD)], @args);
  if (!defined($coverage_threshold)) {
     $coverage_threshold = 0;
  }

  $self->{'coverage_threshold'} = $coverage_threshold;
  my $report = new Bio::Tools::BPlite(@args);

  #start processing
  my $sbjct = $report->nextSbjct;

  if (!defined($sbjct)) {
     return undef;
  }

  $self->{'sbjct_name'}    = $sbjct->name;

  my $sbjct_genomedbrange  = new genomedbrange(-dbrange=>$sbjct->name);
     #my @fields = split(/\s+/, $sbjct->name);
     #my $sbjct_genomedbrange  = new genomedbrange(-dbrange=>$fields[2]);

  if (!defined($sbjct_genomedbrange)) {
     die "sbjct_genomedbrange is not an object of genomedbrange\n";
  }
  $self->{'genomedbrange'} = $sbjct_genomedbrange;

  my $hsp_selected;
  my $hsp_selected_idx;
  my $hsp_idx = -1;

  my @hsps=();
  while (my $hsp = $sbjct->nextHSP) {
     my $coverage = $hsp->query->frac_identical / $hsp->query->seqlength;
     if ($coverage < $coverage_threshold) {
        next;
     }

     $hsp_idx++;
     my $tempStrand = $hsp->query->strand * $hsp->hit->strand;
     
     #seqlevel dbversion seqregion pos1 pos2 strand"
     $hsp->hit->{'genomedbrange'} = new genomedbrange(
                                           -seqlevel =>$sbjct_genomedbrange->getSeqLevel,
                                           -dbversion=>$sbjct_genomedbrange->getDBVersion,
                                           -seqregion=>$sbjct_genomedbrange->getSeqRegion,
                                           -pos1     =>$sbjct_genomedbrange->getPos1 + $hsp->hit->start -1, 
                                           -pos2     =>$sbjct_genomedbrange->getPos1 + $hsp->hit->end -1,
                                           -strand   =>$tempStrand
                                           );

     $hsp->hit->{'coverage'}  = sprintf("%.3f", $coverage);
     $hsp->hit->{'sbjct_dbr'} = $sbjct_genomedbrange;

     if ($tempStrand > 0) {
        my $temp_u_5       = $hsp->query->start - 1;
        my $temp_u_3       = $hsp->query->seqlength - $hsp->query->end;

        my $unaligned_len  = $hsp->query->seqlength - ($hsp->query->end - $hsp->query->start + 1);
        my $suffix_len     = sprintf("%.0f", ($unaligned_len - $temp_u_5 - $temp_u_3 + 1)/2);

        $hsp->hit->{'u_5'} = $temp_u_5 + $suffix_len;
        $hsp->hit->{'u_3'} = $temp_u_3 + $suffix_len;
     } else {
        my $temp_u_3       = $hsp->query->start - 1;
        my $temp_u_5       = $hsp->query->seqlength - $hsp->query->end;
        my $unaligned_len  = $hsp->query->seqlength - ($hsp->query->end - $hsp->query->start + 1);
        my $suffix_len     = sprintf("%.0f", ($unaligned_len - $temp_u_5 - $temp_u_3 + 1)/2);

        $hsp->hit->{'u_3'} = $temp_u_3 + $suffix_len;
        $hsp->hit->{'u_5'} = $temp_u_5 + $suffix_len;
     }

     push(@hsps, $hsp);
     if (defined($hsp_selected)) {
        if ($hsp->bits > $hsp_selected->bits) {
           $hsp_selected     = $hsp;
           $hsp_selected_idx = $hsp_idx;
        }
     } else {
        $hsp_selected     = $hsp;
        $hsp_selected_idx = $hsp_idx;
     }
  }
  
  if ($hsp_idx == -1) {
     return undef;
  }

  @hsps = sort {$b->bits <=> $a->bits} @hsps;

  $self->{'hspcount'}         = $hsp_idx + 1;
  $self->{'hsp_selected'}     = $hsp_selected;
  $self->{'hsps'}             = \@hsps;

  return $self;
}

sub toString {
  my $self = shift;
  return join("\t", $self->numHSPs,
                    $self->bestCoverage,
                    $self->getBestHitGenomeDBRange->getSequenceName,
                    $self->getBestHitExpectedGenomeDBRange->getSequenceName);
}

sub toString2 {
  my $self = shift;
  
  my $tempBestCoverage = $self->bestCoverage;
  my $compatibleCount  = 0;
  
  for (my $i = 0; $i < $self->numHSPs; $i++) {
      if ($self->getHitCoverage($i) >= $tempBestCoverage) {
         #print join("\t", $self->getHitCoverage($i), $tempBestCoverage), "\n";
         $compatibleCount++;
      }
  }
  return $compatibleCount;
}

sub hspSelected {
  my $self = shift;
  return $self->{'hsp_selected'};
}

sub numHSPs {
  my $self = shift;
  return $self->{'hspcount'};
}

sub bestCoverage {
  my $self = shift;
  my $frac_identical = $self->{'hsp_selected'}->query->frac_identical;
  my $query_length   = $self->{'hsp_selected'}->query->seqlength;
  my $coverage = sprintf("%.3f", $frac_identical / $query_length);
  return $coverage;
}

sub getBestHitGenomeDBRange {
  my $self = shift;
  return $self->hspSelected->hit->{'genomedbrange'};
}

sub getBestHitExpectedGenomeDBRange {
  my $self = shift;
  my $bestDBRange = $self->getBestHitGenomeDBRange;
  my $new_range = new genomedbrange(-dbrange=>$bestDBRange->getSequenceName);
  my $bestHSP = $self->hspSelected;

  #print join("\t", $bestHSP->hit->{'u_5'}, $bestHSP->hit->{'u_3'}), "\n";

  $new_range->extend5end($bestHSP->hit->{'u_5'});
  $new_range->extend3end($bestHSP->hit->{'u_3'});

  return $new_range;
}

sub getHitGenomeDBRangeByIdx {
  my $self = shift;
  my $idx  = shift;

  if (!defined($idx)) {$idx = 0;}
  if (!$self->IsValidIdx($idx)) {
     return undef;
  }
  
  return $self->{'hsps'}->[$idx];
}

sub getHitExpectedGenomeDBRange {
  my $self = shift;
  my $idx  = shift;

  if (!defined($idx)) {$idx = 0;}
  if (!$self->IsValidIdx($idx)) {
     return undef;
  }

  my $thisHSP = $self->{'hsps'}->[$idx];
  my $thisGDBR= $thisHSP->hit->{'genomedbrange'};
  
  my $new_range = new genomedbrange(-dbrange=>$thisGDBR->getSequenceName);
  $new_range->extend5end($thisHSP->hit->{'u_5'});
  $new_range->extend3end($thisHSP->hit->{'u_3'});
  
  return $new_range;
}

sub getHitCoverage {
  my $self = shift;
  my $idx  = shift; 
  
  if (!defined($idx)) {$idx = 0;}
  if (!$self->IsValidIdx($idx)) {
     return undef;
  }

  return $self->{'hsps'}->[$idx]->hit->{'coverage'};
}

sub IsValidIdx {
  my $self = shift;
  my $idx  = shift;
  
  if ($idx > $self->numHSPs - 1) {
     warn "idx outofarraybound\n";
     return 0;
  } else {
     return 1;
  }
}

1;
