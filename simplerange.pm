package simplerange;
use math;
use strict;
use Bio::Root::Root;
use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'simplerange';
$VERSION  = 0.1;

=head1 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($pos1, $pos2, $seqname, $strand, $desc, $geneclass) = 
     $self->_rearrange( [qw(POS1 POS2 SEQNAME STRAND DESC GENECLASS)], @args);
  
  if ($pos1 > $pos2) {
     if (defined($strand)) {
        warn "$pos1 is greater than $pos2, however, you've defined strand as $strand\n";
        return undef;
     } else {
        $strand = -1;
        my $temp = $pos1; 
        $pos1 = $pos2;
        $pos2 = $temp;
     }
  }
  
  if (!defined($strand)) {
     $strand = 1;
  }

  $self->{"pos1"}       = $pos1;
  $self->{"pos2"}       = $pos2;
  $self->{"seqname"}    = $seqname;
  $self->{"strand"}     = $strand;
  $self->{"geneclass"}  = $geneclass;
  $self->{"desc"}       = $desc;

  return $self; # success - we hope!
}

sub setPos1 {
  my $self    = shift;
  my $newPos1 = shift;
  $self->{'pos1'} = $newPos1;
}

sub setPos2 {
  my $self = shift;
  my $newPos2 = shift;
  $self->{'pos2'} = $newPos2;
}

sub extend5end {
  my $self = shift;
  my $distance = shift;
  if ($distance < 0) {
     die "Distance $distance is less than zero\n";
  }
  
  $self->setPos1($self->getPos1 - $distance);
}

sub extend3end {
  my $self = shift;
  my $distance = shift;
  if ($distance < 0) {
     die "Distance $distance is less than zero\n";
  }

  $self->setPos2($self->getPos2 + $distance);
}


sub getSeqName {
  my $self = shift;
  return $self->{'seqname'};
}

sub getLen {
  my $self = shift;
  my $len  = $self->getPos2 - $self->getPos1 + 1;
  return $len;
}

sub getDownStreamLen {
  my $self = shift;
  my $len  = shift;

  if (!defined($len) || $len <=0) {
     return undef;
  }

  my $new_pos1;
  my $new_pos2;
  if ($self->getStrand == -1) {
     $new_pos1 = $self->getPos1 - $len;
     $new_pos2 = $self->getPos1 - 1;
  } else {
     $new_pos1 = $self->getPos2 + 1;
     $new_pos2 = $self->getPos2 + $len;
  }

  return new simplerange(-pos1=>$new_pos1, 
                         -pos2=>$new_pos2, 
                         -seqname=>$self->getSeqName, 
                         -strand=>$self->getStrand,
                         -desc  =>$self->getDesc
                        );
}


sub getUpStreamLen {
  my $self = shift;
  my $len  = shift;

  if (!defined($len) || $len <=0) {
     return undef;
  }

  my $new_pos1;
  my $new_pos2;
  if ($self->getStrand == -1) {
     $new_pos1 = $self->getPos2 + 1;
     $new_pos2 = $self->getPos2 + $len;
  } else {
     $new_pos1 = $self->getPos1 - $len ;
     $new_pos2 = $self->getPos2 - 1;
  }

  return new simplerange(-pos1=>$new_pos1, 
                         -pos2=>$new_pos2, 
                         -seqname=>$self->getSeqName, 
                         -strand=>$self->getStrand,
                         -desc  =>$self->getDesc
                        );
}

sub getPos1 {
  my $self = shift;
  return $self->{"pos1"};
}

sub getPos2 {
  my $self = shift;
  return $self->{"pos2"};
}

sub getgeneclass {
  my $self = shift;
  return $self->{'geneclass'};
}

sub isContained {
  my $self = shift;
  my $tempRange = shift;

  if (!$tempRange->isa('simplerange')) {
     warn "The range you've passed is not an object of simplerange\n";
     return undef;
  }

  my $overlap_start = max($tempRange->getPos1, $self->getPos1);
  my $overlap_end   = min($tempRange->getPos2, $self->getPos2);
  my $overlap_len   = $overlap_end - $overlap_start + 1;
  
  if ($overlap_len == $self->getLen || $overlap_len == $tempRange->getLen) {
     return 1;
  } else {
     return undef;
  }
}

sub checkOverlap {
  my $self = shift;
  my $temprange = shift;
  
  if (!$temprange->isa("simplerange")) {
     warn "The range you've passed is not an object of simplerange\n";
     return undef;
  }

  my $overlap_start = max($temprange->getPos1, $self->getPos1);
  my $overlap_end   = min($temprange->getPos2, $self->getPos2);
  my $overlap_len   = $overlap_end - $overlap_start + 1;
 
  if ($overlap_len > 0) {
     return new simplerange(-pos1=>$overlap_start, -pos2=>$overlap_end, -seqname=>"overlap_check");
  } else {
     return undef;
  }
}

sub checkOverlapTargetExt3 {
  my $self      = shift;
  my $temprange = shift;
  my $ext3      = shift;
  
  if (!$temprange->isa("simplerange")) {
     warn "The range you've passed is not an object of simplerange\n";
     return undef;
  }

  my $overlap_start = max($temprange->getPos1, $self->getPos1);
  my $overlap_end   = min($temprange->getPos2 + $ext3, $self->getPos2);
  my $overlap_len   = $overlap_end - $overlap_start + 1;
 
  if ($overlap_len > 0) {
     return new simplerange(-pos1=>$overlap_start, -pos2=>$overlap_end, -seqname=>"overlap_check");
  } else {
     return undef;
  }
}

sub checkOverlapTargetExt5 {
  my $self      = shift;
  my $temprange = shift;
  my $ext5      = shift;
  
  if (!$temprange->isa("simplerange")) {
     warn "The range you've passed is not an object of simplerange\n";
     return undef;
  }

  my $overlap_start = max($temprange->getPos1 - $ext5, $self->getPos1);
  my $overlap_end   = min($temprange->getPos2, $self->getPos2);
  my $overlap_len   = $overlap_end - $overlap_start + 1;

  #print join("\t",$overlap_start, $overlap_end, $overlap_len), "\n";
 
  if ($overlap_len > 0) {
     return new simplerange(-pos1=>$overlap_start, -pos2=>$overlap_end, -seqname=>"overlap_check");
  } else {
     return undef;
  }
}

sub checkOverlapRaw {
  my $self = shift;
  my $pos1 = shift;
  my $pos2 = shift;
  
  my $new_range = new simplerange(-POS1=>$pos1, -POS2=>$pos2);
  return $self->checkOverlap($new_range);
}

sub mergeRange {
  my $self = shift;
  my $temprange = shift;
  
  if (!$temprange->isa("simplerange")) {
     warn "The range you've passed is not an object of simplerange\n";
     return undef;
  }

  my $new_start  = min($temprange->getPos1, $self->getPos1);
  my $new_end    = max($temprange->getPos2, $self->getPos2);
  my $new_strand = $self->getStrand;
 
  if ($temprange->getStrand != $new_strand) {
     $new_strand = 0;
  }

  return new simplerange(-pos1=>$new_start, -pos2=>$new_end, -strand=>$new_strand, -seqname=>"merged");
}

sub mergeOverlappedRange {
  my $self = shift;
  my $temprange = shift;
  
  if (defined($self->checkOverlap($temprange))) {
     return $self->mergeRange($temprange);
  } else {
     return undef;
  }  
}

sub getSubRange {
  my $self = shift; 
  my $subpos1 = shift;
  my $subpos2 = shift;

  if ($subpos1 >= $subpos2) {
     return undef;
  }
  
  my $newR;

  if ($self->getStrand != -1) {
     $newR = new simplerange(-seqname => $self->getSeqName, 
                             -pos1    => $self->getPos1 + $subpos1 -1,
                             -pos2    => $self->getPos1 + $subpos2 -1,
                             -strand  => $self->getStrand
                            );
  } else {
     $newR = new simplerange(-seqname =>$self->getSeqName,
                             -pos1    =>$self->getPos2 - $subpos2 + 1,
                             -pos2    =>$self->getPos2 - $subpos1 + 1,
                             -strand  =>$self->getStrand
                            );
  }
  
  return $newR;
}

sub getDesc {
  my $self = shift;
  return $self->{"desc"};
}

sub getStrand {
  my $self = shift;
  return $self->{"strand"};
}

sub toString {
  my $self = shift;
  return join(":", $self->getPos1, $self->getPos2, $self->getStrand);
}

sub toString2 {
  my $self = shift;
  return join(":", $self->getgeneclass, $self->getPos1, $self->getPos2, $self->getStrand);
}

sub toString3 {
  my $self = shift;
  return join(":", $self->getSeqName, $self->getPos1, $self->getPos2, $self->getStrand);
}

sub compare {
  my ($self, $range2) = (shift, shift);
  if (!$range2->isa("simplerange")) {die "The input range is not a simplerange\n";}
  return ($self->getPos1 <=> $range2->getPos1 ||
          $self->getPos2 <=> $range2->getPos2);
}

sub compare1 {
  my ($self, $range2) = (shift, shift);
  if (!$range2->isa("simplerange")) {die "The input range is not a simplerange\n";}
  return ($self->getPos1 <=> $range2->getPos1);
  
}

1;
