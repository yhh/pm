package genomedbrange;
use strict;
use Bio::Root::Root;
use chr_format;
use math;
use string_hash_code;
use simplerange;
use sortfuncs;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'genomedbrange';
$VERSION  = 0.1;


=head1 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($sequence_name, $seqLevel, $dbVersion, $seqRegion, $pos1, $pos2, $strand) = 
      $self->_rearrange([qw(DBRANGE SEQLEVEL DBVERSION SEQREGION POS1 POS2 STRAND)], @args);
  
  if (defined($sequence_name)) {
     $self->{'sequence_name'} = $sequence_name;
     if ($sequence_name!~/\:/) {
        return undef;
     }
     $sequence_name=~s/ //g;
     ($seqLevel, $dbVersion, $seqRegion, $pos1, $pos2, $strand) = split(":", $sequence_name);
  } else {
     if (!defined($seqRegion)) {
        return undef;
     }
     
     if (!defined($dbVersion)) {
        $dbVersion = 'NCBI35';
     }
     
     if (!defined($seqLevel)) {
        $seqLevel = 'chromosome';
     }

     if ($pos1 > $pos2) {
        if (defined($strand)) {
           warn "$pos1 is greater than $pos2, however, you've defined strand $strand\n";
           return undef;
        }
        my $temp = $pos1; 
        $pos1 = $pos2;
        $pos2 = $temp;
        $strand = -1;
     } elsif (!defined($strand)) {
        $strand = 1;
     }
     $sequence_name = join(':', $seqLevel, $dbVersion, $seqRegion, $pos1, $pos2, $strand);
     $sequence_name=~s/\s//g;
     $self->{'sequence_name'} = $sequence_name;
  }

  my ($sub_pos1, $sub_pos2) = (0, 0);
  if ($sequence_name =~/\//) {
     my ($temp1, $temp2) = split("\/", $sequence_name);
     $sequence_name = $temp1;
     ($sub_pos1, $sub_pos2) = split("-", $temp2);
  }

  $self->{"seqLevel"}  = $seqLevel;
  $self->{"dbVersion"} = $dbVersion;
  $self->{"seqRegion"} = $seqRegion;
  if ($sub_pos1 > 0) {
     $self->{"pos1"}   = $pos1 + $sub_pos1 -1;
  } else {
     $self->{"pos1"}   = $pos1;
  }
  if ($sub_pos2 > 0) {
     $self->{"pos2"}   = $pos1 + $sub_pos2 -1;
  } else {
     $self->{"pos2"}   = $pos2;
  }
  $self->{"strand"}    = $strand;
  
  #chromosome:NCBIM32:5:29207146:29225145:1/1501-1800
  return $self; # success - we hope!
}

sub getBoundDBR {
  my $self  = shift;
  my $bGDBR = shift;

  if (!$bGDBR->isa("genomedbrange")) {
     return undef;
  }

  if (isCompareCompatible($self, $bGDBR)) {
     
     my @ordered_GDBRs = sort sortfuncs::orderbygenomedbrange ($self, $bGDBR);

     #print join("\t", $ordered_GDBRs[0]->getSequenceName, 
     #                 $ordered_GDBRs[1]->getSequenceName
     #          ), "\n";

     my $newPos1 = $ordered_GDBRs[0]->getPos2 + 1;
     my $newPos2 = $ordered_GDBRs[1]->getPos1 - 1;

     if ($newPos1 > $newPos2) {
        return undef;
     }

     my $new_gdbr = new genomedbrange(-SEQLEVEL   => $self->getSeqLevel, 
                                      -DBVERSION  => $self->getDBVersion,
                                      -SEQREGION  => $self->getSeqRegion,
                                      -POS1       => $newPos1,
                                      -POS2       => $newPos2
                                     );
     return $new_gdbr;
  } else {
     return undef;
  }
}

sub mergeDBRStrand {
  my $self  = shift;
  my $bGDBR = shift;

  if (!$bGDBR->isa("genomedbrange")) {
     return undef;
  }

  if (isCompareCompatible($self, $bGDBR)) {
     
     my $newPos1   = min($self->getPos1, $bGDBR->getPos1);
     my $newPos2   = max($self->getPos2, $bGDBR->getPos2);
     my $newStrand = $self->getStrand;

     my $new_gdbr = new genomedbrange(-SEQLEVEL   => $self->getSeqLevel, 
                                      -DBVERSION  => $self->getDBVersion,
                                      -SEQREGION  => $self->getSeqRegion,
                                      -POS1       => $newPos1,
                                      -POS2       => $newPos2,
                                      -STRAND 	  => $newStrand
                                     );
     return $new_gdbr;
  } else {
     return undef;
  }	
}

sub mergeDBR {
  my $self  = shift;
  my $bGDBR = shift;

  if (!$bGDBR->isa("genomedbrange")) {
     return undef;
  }

  if (isCompareCompatible($self, $bGDBR)) {
     
     my $newPos1 = min($self->getPos1, $bGDBR->getPos1);
     my $newPos2 = max($self->getPos2, $bGDBR->getPos2);

     my $new_gdbr = new genomedbrange(-SEQLEVEL   => $self->getSeqLevel, 
                                      -DBVERSION  => $self->getDBVersion,
                                      -SEQREGION  => $self->getSeqRegion,
                                      -POS1       => $newPos1,
                                      -POS2       => $newPos2
                                     );
     return $new_gdbr;
  } else {
     return undef;
  }
}

sub isCompareCompatible {
  my $aGDBR = shift;
  my $bGDBR = shift;
  
  if ($aGDBR->getSeqLevel  ne $bGDBR->getSeqLevel  || 
      $aGDBR->getDBVersion ne $bGDBR->getDBVersion ||
      $aGDBR->getSeqRegion ne $bGDBR->getSeqRegion) {
      return 0;
  }
  return 1;
}

sub getDownStreamLen {
  my $self = shift;
  my $len  = shift;
  if (!defined($len) || $len <= 0 ) {
     return undef;
  } else {
     my $new_pos1;
     my $new_pos2;
     if ($self->getStrand == -1) {
        $new_pos1 = $self->getPos1 - $len;
        $new_pos2 = $self->getPos1 - 1;
     } else {
        $new_pos1 = $self->getPos2 + 1;
        $new_pos2 = $self->getPos2 + $len
     }

     return new genomedbrange(-seqlevel =>$self->getSeqLevel,
                              -dbversion=>$self->getDBVersion,
                              -seqregion=>$self->getSeqRegion,
                              -pos1     =>$new_pos1, 
                              -pos2     =>$new_pos2,
                              -strand   =>$self->getStrand);
  }
}

sub hashCode {
  my $self = shift;
  return join('.', $self->getSeqRegion, $self->getPos1, $self->getPos2, $self->getStrand);
}

sub getDownStreamRange {
  my $self    = shift;
  my $range1  = shift;
  my $range2  = shift;

  my $len     = $range2 - $range1 + 1;

  if (!defined($range1) || !defined($range2) || $len < 1) {
     return undef;
  } else {
     my $new_pos1;
     my $new_pos2;
     if ($self->getStrand == -1) {
        $new_pos1 = $self->getPos1 - $range2 + 1;
        $new_pos2 = $self->getPos1 - $range1;
     } else {
        $new_pos1 = $self->getPos2 + $range1;
        $new_pos2 = $self->getPos2 + $range2 - 1;
     }

     return new genomedbrange(-seqlevel =>$self->getSeqLevel,
                              -dbversion=>$self->getDBVersion,
                              -seqregion=>$self->getSeqRegion,
                              -pos1     =>$new_pos1, 
                              -pos2     =>$new_pos2,
                              -strand   =>$self->getStrand);
  }
}

sub getLength {
  my $self = shift;
  return ($self->getPos2 - $self->getPos1 + 1);
}

sub relativeTo {
  my $self = shift;
  my $dbr2 = shift;
  if (!$dbr2->isa("genomedbrange")) {
     my $new_dbr2 = new genomedbrange(-dbrange=>$dbr2);
     if (!defined($new_dbr2)) {
        die "Cannot initiate $dbr2\n";
     }
     $dbr2 = $new_dbr2;
  }

  if ($dbr2->getStrand != $self->getStrand) {
     warn "strand\n";
     if ($dbr2->getStrand == -1) {
        return new simplerange(-pos1=>($dbr2->getPos1 - $self->getPos1) + 1,
                               -pos2=>($dbr2->getPos2 - $self->getPos1) + 1,
                               -strand=>1,
                               -seqname=>$self->getSequenceName
                              );
     } else {
        return undef;
     }
  } else {
     if ($dbr2->getStrand == -1) {
         return new simplerange(-pos1=>abs($dbr2->getPos2 - $self->getPos2) + 1, 
                                -pos2=>abs($dbr2->getPos1 - $self->getPos2) + 1, 
                                -strand=>1,
                                -seqname=>$self->getSequenceName);
     } else {
         return new simplerange(-pos1=>   ($dbr2->getPos1 - $self->getPos1) + 1, 
                                -pos2=>   ($dbr2->getPos2 - $self->getPos1) + 1, 
                                -strand=>1,
                                -seqname=>$self->getSequenceName);
     }
  }
}

sub getUpStreamLen {
  my $self = shift;
  my $len  = shift;
  if (!defined($len) || $len <= 0 ) {
     return undef;
  } else {
     my $new_pos1;
     my $new_pos2;
     if ($self->getStrand == -1) {
        $new_pos1 = $self->getPos2 + 1;
        $new_pos2 = $self->getPos2 + $len;
     } else {
        $new_pos1 = $self->getPos1 - $len;
        $new_pos2 = $self->getPos1 - 1;
     }

     return new genomedbrange(-seqlevel =>$self->getSeqLevel,
                              -dbversion=>$self->getDBVersion,
                              -seqregion=>$self->getSeqRegion,
                              -pos1     =>$new_pos1, 
                              -pos2     =>$new_pos2,
                              -strand   =>$self->getStrand);
  }
}

sub getSequenceName {
  my $self = shift;
  return $self->{'sequence_name'};
}

sub getSequenceName2 {
  my $self = shift;
  my $name = $self->getSequenceName;
  $name=~s/\:/\./g;
  return $name;
}

sub compare {
  my ($self, $range2) = (shift, shift);
  if (!$range2->isa("genomedbrange")) {die "The input range is not a genomedbrange\n";}
  if ( $range2->getDBVersion ne $self->getDBVersion) {die "The DB versions are not compatible\n";}
  return (hashcode($self->getSeqRegion) <=> hashcode($range2->getSeqRegion) ||
          $self->getPos1 <=> $range2->getPos1 ||
          $self->getPos2 <=> $range2->getPos2);
}

sub getSimpleRange {
  my $self = shift;
  my $res_simple_range = new simplerange(-pos1=>$self->getPos1,  -pos2=>$self->getPos2, 
                                         -seqname=>$self->getSeqRegion, -strand=>$self->getStrand);
  return $res_simple_range;
}

sub checkOverlap {
  my ($self, $range2) = (shift, shift);
  if (!$range2->isa("genomedbrange")) {die "The input range is not a genomedbrange\n";}
  
  if ($self->getSeqRegion() ne $range2->getSeqRegion()) {
     return undef;
  }

  if ($range2->getDBVersion ne $self->getDBVersion) {die "The DB versions are not compatible\n";}
  
  my $overlap_start = max($self->getPos1, $range2->getPos1);
  my $overlap_end   = min($self->getPos2, $range2->getPos2);
  my $overlap_length= $overlap_end - $overlap_start + 1;
  
  if ($overlap_length > 0) {
     my $new_range = new genomedbrange(-dbrange=>$self->getSequenceName);

     $new_range->setPos1($overlap_start);
     $new_range->setPos2($overlap_end);
     
     if ($self->getStrand == -1) {
        warn "All overlap gdbr are strand free\n";
        $new_range->reverseStrand;
     }

     return $new_range;
  } else {
     return undef;
  }
}

sub checkOverlapRaw {
  my ($self, $rawrange2) = (shift, shift);
  my $range2 = new genomedbrange(-dbrange=>$rawrange2);
  if (!defined($range2)) {
     die "Please check your range: $rawrange2\n";
  }

  return $self->checkOverlap($range2);
}

sub compareRaw {
  my ($self, $rawRange2) = (shift, shift);
  my $range2 = new genomedbrange(-dbrange=>$rawRange2);
  if (!$range2->isa("genomedbrange")) {die "The input range is not a genomedbrange\n";}
  if ( $range2->getDBVersion ne $self->getDBVersion) {
     warn join("\t", $self->getSequenceName, $range2->getSequenceName), "\n";
     die "The DB versions are not compatible:".$range2->getDBVersion." ".$self->getDBVersion."\n";
  }
  return (hashcode($self->getSeqRegion) <=> hashcode($range2->getSeqRegion) ||
          $self->getPos1 <=> $range2->getPos1 ||
          $self->getPos2 <=> $range2->getPos2);
}

sub getDistance {
  my $self   = shift;
  my $range2 = shift;
  if (!$range2->isa("genomedbrange")) {die "The input range is not a genomedbrange\n";}
  if ( $range2->getDBVersion ne $self->getDBVersion) {
     warn join("\t", $self->getSequenceName, $range2->getSequenceName), "\n";
     die "The DB versions are not compatible\n";
  }
  
  if ( $range2->getSeqRegion ne $self->getSeqRegion) {return undef;}
  return abs($range2->getPos1 - $self->getPos1);
}

sub getDistance2 {
  my $self   = shift;
  my $range2 = shift;

  if (!$range2->isa("genomedbrange")) {die "The input range is not a genomedbrange\n";}
  if ( $range2->getDBVersion ne $self->getDBVersion) {
     warn join("\t", $self->getSequenceName, $range2->getSequenceName), "\n";
     die "The DB versions are not compatible:".$range2->getDBVersion." ".$self->getDBVersion."\n";
  }

  if ($range2->getSeqRegion ne $self->getSeqRegion) {
     return undef;
  }
  
  if ($range2->getPos2 <= $self->getPos1) {
     return abs($range2->getPos2 - $self->getPos1);
  } elsif ($range2->getPos1 >= $self->getPos2) {
     return abs($range2->getPos1 - $self->getPos2);
  } else {
     return 0;
  }
}

sub clone {
  my $self = shift;
  return new genomedbrange(-SEQLEVEL   => $self->getSeqLevel,
                           -DBVERSION  => $self->getDBVersion,
                           -SEQREGION  => $self->getSeqRegion,
                           -POS1       => $self->getPos1,
                           -POS2       => $self->getPos2,
                           -STRAND     => $self->getStrand
                          );
}

sub shiftDistance {
  my $self = shift;
  my $distance = shift;
  my $direction = shift;

  my $new_range = new genomedbrange(-dbrange=>$self->getSequenceName);
  my $shifted_distance = $distance * $direction; 

  $new_range->setPos1($self->getPos1 + $shifted_distance);
  $new_range->setPos2($self->getPos2 + $shifted_distance);

  return $new_range;
}

sub getSubRange {
  my $self = shift;
  my $pos1 = shift;
  my $pos2 = shift;
  my $strand = shift;
  
  if (!defined($pos1)) {
     $pos1 = 1;
  }

  if (!defined($pos2)) {
     $pos2 = $self->getLength;
  }

  if (!defined($strand)) {
     $strand = 1;
  }

  my $new_range = new genomedbrange(-dbrange=>$self->getSequenceName);

  if ($self->getStrand == -1) {
     #print join("\t", $pos1, $pos2), "\n";
     $new_range->setPos1($self->getPos2 - $pos2 + 1);
     $new_range->setPos2($self->getPos2 - $pos1 + 1);
  } else {
     #print join("\t", $pos1, $pos2), "\n";
     $new_range->setPos1($self->getPos1 + $pos1 - 1);
     $new_range->setPos2($self->getPos1 + $pos2 - 1);
  }

  if ($strand == -1) {
    $new_range->reverseStrand;
  }

  return $new_range;
}

sub extend5end {
  my $self = shift;
  my $distance = shift;
  if ($distance < 0) {
     die "Distance $distance is less than zero\n";
  }
  
  $self->setPos1($self->getPos1 - $distance);
}

sub extend5end_directional {
  my $self = shift;
  my $distance = shift;
  if ($distance < 0) {
     die "Distance $distance is less than zero\n";
  }
  
  if ($self->getStrand < 0) {
     $self->setPos2($self->getPos2 + $distance);
  } else {
     $self->setPos1($self->getPos1 - $distance);
  }
}

sub extend3end {
  my $self = shift;
  my $distance  = shift;

  if ($distance < 0) {
     die "Distance $distance is less than zero\n";
  }

  $self->setPos2($self->getPos2 + $distance);
}

sub extend3end_directional {
  my $self = shift;
  my $distance  = shift;

  if ($distance < 0) {
     die "Distance $distance is less than zero\n";
  }

  if ($self->getStrand < 0) {
     $self->setPos1($self->getPos1 - $distance);
  } else {
     $self->setPos2($self->getPos2 + $distance);
  }
}

sub getDBVersion {
  my $self = shift;
  return $self->{"dbVersion"};
}

sub getSeqLevel {
  my $self = shift;
  return $self->{"seqLevel"};
}

sub getSeqRegion {
  my $self = shift;
  return $self->{"seqRegion"};
}

sub getPos1 {
  my $self = shift;
  return $self->{"pos1"};
}

sub setPos1 {
  my $self = shift;
  my $newPos1 = shift;
  $self->{'pos1'} = $newPos1;
  my $sequence_name = join(':', $self->getSeqLevel, 
                                $self->getDBVersion, 
                                $self->getSeqRegion, 
                                $self->getPos1, 
                                $self->getPos2, 
                                $self->getStrand);
  
  $self->{'sequence_name'} = $sequence_name;
}

sub getPos2 {
  my $self = shift;
  return $self->{"pos2"};
}

sub setPos2 {
  my $self = shift;
  my $newPos2 = shift;
  $self->{'pos2'} = $newPos2;
  my $sequence_name = join(':', $self->getSeqLevel, 
                                $self->getDBVersion, 
                                $self->getSeqRegion, 
                                $self->getPos1,
                                $self->getPos2, 
                                $self->getStrand);
  $self->{'sequence_name'} = $sequence_name;
}

sub getStrand {
  my $self = shift;
  return $self->{"strand"};
}

sub setStrand {
	my $self = shift;
	my $strand = shift;
	if (!defined($strand)) {
		$strand = 1;
	}
	$self->{'strand'} = $strand;
}

sub reverseStrand {
  my $self = shift;
  $self->{'strand'} = $self->getStrand * -1;
  my $sequence_name = join(':', $self->getSeqLevel,
                                $self->getDBVersion,
                                $self->getSeqRegion,
                                $self->getPos1,
                                $self->getPos2,
                                $self->getStrand);
  $self->{'sequence_name'} = $sequence_name;
}

sub get_chr_format {
  my $self = shift;
  #CHR_NAME POS1 POS2 STRAND
  return chr_format->newChrPos(
                      -chr_name=>$self->getSeqRegion,
                      -pos1    =>$self->getPos1,
                      -pos2    =>$self->getPos2,
                      -strand  =>$self->getStrand
                   );
}

1;
