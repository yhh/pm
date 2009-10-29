package GDBRList;

use strict;
use Bio::Root::Root;
use genomedbrange;
use sortfuncs;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'GDBRList';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  return $self; # success - we hope!
}

sub add {
  my $self = shift;
  my $gdbr = shift;
  my $obj  = shift;
  
  if (!defined($self->{'hash'}->{$gdbr->getSequenceName})) {
     $self->{'hash'}->{$gdbr->getSequenceName} = $obj;
  } else {
     warn "cannot add object because ".$gdbr->getSequenceName." already exists\n";
  }

  $self->{'orderedRawGDBRs'} = undef;
}

sub getObjByRawGDBR {
  my $self = shift;
  my $GDBR_raw = shift;
  
  return $self->{'hash'}->{$GDBR_raw};
}

sub getObjByGDBR {
  my $self = shift;
  my $GDBR = shift;
  
  return $self->{'hash'}->{$GDBR->getSequenceName};
}

sub getObjAt {
  my $self = shift;
  my $idx  = shift;

  my $rawGDBR = undef;

  if (defined($self->{'orderedRawGDBRs'})) {
     $rawGDBR = $self->{'orderedRawGDBRs'}->[$idx];
  }

  if (defined($rawGDBR)) {
     return $self->getObjByRawGDBR($rawGDBR);
  } else {
     return undef;
  }
}

sub getOrderedGDBR {
  my $self = shift;
  my @keys = keys %{$self->{'hash'}};

  my @ordered_GDBRs = sort sortfuncs::orderbygenomedbrange_raw @keys ;

  $self->{'orderedRawGDBRs'} = \@ordered_GDBRs;
  return \@ordered_GDBRs;
}

sub findRange {
  my $self = shift;
  my $dbr  = shift;

  my $lowerIdx  = $self->dividefinder($dbr);
  #print $lowerIdx, "\n";
  my $lowerGDBR = new genomedbrange(-dbrange=>$self->{'orderedRawGDBRs'}->[$lowerIdx]);

  #print $lowerGDBR->getSequenceName, "kkkk\n";

  if ($dbr->getSeqRegion ne $lowerGDBR->getSeqRegion) {
     return undef, undef;
  }

  my $newPUpperIdx = $lowerIdx;

  while ($newPUpperIdx >= 0 && $newPUpperIdx < scalar(@{$self->{'orderedRawGDBRs'}})) {
     my $upperRawGDBR = $self->{'orderedRawGDBRs'}->[$newPUpperIdx];
     my $upperGDBR = new genomedbrange(-dbrange=>$upperRawGDBR);
     
     if (!defined($dbr->checkOverlapRaw($upperRawGDBR))) {
        $newPUpperIdx--;
        last;
     }
     $newPUpperIdx++;
  }

  if ($newPUpperIdx >= scalar(@{$self->{'orderedRawGDBRs'}})) {
     $newPUpperIdx = scalar(@{$self->{'orderedRawGDBRs'}}) - 1;
  }

  if ($newPUpperIdx >= $lowerIdx) {
     return $lowerIdx, $newPUpperIdx;
  } else {
     return undef, undef;
  }
}

sub dividefinder {
  my $self = shift;
  my ($dbr, $lower, $upper) = (shift, shift, shift);

  if (!defined($lower)) {
     $lower = 0;
  }

  if (!defined($self->{'orderedRawGDBRs'})) {
     $self->getOrderedGDBR();
  }

  if (!defined($upper)) {
     $upper = scalar(@{$self->{'orderedRawGDBRs'}})-1;
  }

  if ($upper - $lower == 1) {
     if (defined($dbr->checkOverlapRaw($self->{'orderedRawGDBRs'}->[$lower]))) {
        return $lower;
     } else {
        return $upper;
     }
  }

  my $middle=int(($upper+$lower)/2);
  #print join("\t", $upper, $middle, $lower), " upper middle lower\n";
  #print $self->{'orderedRawGDBRs'}->[$lower], "\n";
  #if (($dbr->compareRaw($self->{'orderedRawGDBRs'}->[$lower])  >=0) &&
  #    ($dbr->compareRaw($self->{'orderedRawGDBRs'}->[$middle]) < 0)) {
  if ($dbr->compareRaw($self->{'orderedRawGDBRs'}->[$middle]) < 0) {
     return $self->dividefinder($dbr, $lower, $middle);
  } else {
     return $self->dividefinder($dbr, $middle, $upper);
  }
}

1;
