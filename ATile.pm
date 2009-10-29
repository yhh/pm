package ATile;

use strict;
use Bio::Root::Root;
use math;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'ATile';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  ($self->{'start'}, $self->{'end'}, $self->{'strand'}, $self->{'dbrange'}) = 
     $self->_rearrange([qw(START END STRAND DBRANGE)], @args);

  if (!defined($self->{'start'})) {
     $self = undef;
  }

  if (!defined($self->{'end'})) {
     $self = undef;
  }

  if (defined($self)) {
     if ($self->{'start'} > $self->{'end'}) {
        my $temp = $self->{'end'};
        $self->{'end'}  = $self->{'start'};
        $self->{'start'}= $temp;
     }
  }
 
  if (!defined($self->{'strand'})) {
     $self->{'strand'} = 1;
  }

  return $self; # success - we hope!
}

sub start {
  my $self = shift;
  return $self->{'start'};
}

sub end {
  my $self = shift;
  return $self->{'end'};
}

sub length {
  my $self = shift;

  return ($self->end - $self->start + 1);
}

sub strand {
  my $self = shift;
  return $self->{'strand'};
}

sub isOverlap {
  my $self = shift;

  my ($frag_start, $frag_end) = sortNum(shift, shift);

  my $start = $self->{'start'};
  my $end   = $self->{'end'};

  my $overlap_start = max($start, $frag_start);
  my $overlap_end   = min($end,   $frag_end);

  my $overlap_len   = $overlap_end - $overlap_start + 1;

  return ($overlap_len > 0);
}

sub isContained {
  my $self = shift;
  my ($frag_start, $frag_end) = sortNum(shift, shift);

  my $start = $self->{'start'};
  my $end   = $self->{'end'};

  return ($start <= $frag_start && $end >= $frag_end);
}

sub overlapLen {
  my $self = shift;
  my ($frag_start, $frag_end) = sortNum(shift, shift);
  
  my $start = $self->{'start'};
  my $end   = $self->{'end'};

  my $overlap_start = max($start, $frag_start);
  my $overlap_end   = min($end,   $frag_end);

  my $overlap_len   = $overlap_end - $overlap_start + 1;
  
  return $overlap_len;
}

sub subDBRange {
  my $self = shift;
  my ($frag_start, $frag_end) = (shift, shift);

  if ($self->dbRange->getStrand < 0) {
     die "The strand of the original genomedbrange cannot be negative\n";
  }

  my $strand = 1;
  if ($frag_start > $frag_end) {
     my $temp = $frag_start;
     $frag_start = $frag_end;
     $frag_end   = $temp;
     $strand     = -1;
  }
  
  my $dist5;
  my $dist3;

  if ($self->strand < 0) {
     $dist5 = $self->end - $frag_end + 1;
     $dist3 = $self->end - $frag_start + 1;
  } else {
     $dist5 = $frag_start - $self->start + 1;
     $dist3 = $frag_end   - $self->start + 1;
  }

  my $newdbrange = $self->dbRange->getSubRange($dist5, $dist3, $strand * $self->strand);
  return $newdbrange;
}

sub dbRange {
  my $self = shift;
  return $self->{'dbrange'};
}

sub sortNum {
  my $num1 = shift;
  my $num2 = shift;
  
  if ($num1 > $num2) {
     my $temp = $num2;
     $num2    = $num1;
     $num1    = $temp;
  }

  return $num1, $num2;
}


