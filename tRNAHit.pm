package tRNAHit;

use strict;
use Bio::Root::Root;
use simplerange;

use vars qw($ID $VERSION @ISA $DESC);

@ISA  = qw(Bio::Root::Root);
$ID   = 'tRNAHit';
$VERSION  = 0.1;
$DESC = "only for use together with eufindtRNA";

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{'seqname'} = $self->_rearrange([qw(SEQNAME)], @args);
  return $self; # success - we hope!
}

sub getBBox {
  my $self = shift;
  return $self->{'BBox'};
}

sub IsGoodHit {
  my $self = shift;
  if (scalar(@_) > 0) {
     $self->{'IsGood'} = shift @_;
  } else {
     if (!defined($self->{'IsGood'})) {
        $self->{'IsGood'} = 0;
     }
     return $self->{'IsGood'};
  }
}

sub addBBox {
  my $self = shift;
  my $line_for_b_box = shift;
  #Bbox at 48 (End=70), Sc= -2.17
  if ($line_for_b_box =~/Bbox at (\d+).*?Sc= (.*)/) {
     my $pos1  = $1;
     my $score = $2;
     my $pos2  = $pos1 + 10;
     my $aRange = new simplerange(-pos1=>$pos1, -pos2=>$pos2, -desc=>$score);
     $self->{'BBox'} = $aRange;
     return 1;
  } else {
     return 0
  }
}

sub addABox {
  my $self = shift;
  my $line_for_a_box = shift;
  #Abox at 2 (St=-3) A:-61.45 AB(25):-0.46 I:-64.08
  if ($line_for_a_box =~/Abox at (\d+).*?A.(.*?)\s/) {
     my $pos1  = $1;
     my $score = $2;
     my $pos2  = $pos1 + 20;
     my $aRange = new simplerange(-pos1=>$pos1, -pos2=>$pos2, -desc=>$score);
     push(@{$self->{'ABox'}}, $aRange);
     if (defined($self->{'BestABox'})) {
        if ($score > $self->{'BestABox'}->getDesc) {
           $self->{'BestABox'} = $aRange;
        }
     } else {
        $self->{'BestABox'} = $aRange;
     }
     return 1;
  } else {
     return 0;
  }
}

sub getBestABox {
  my $self = shift;
  return $self->{'BestABox'};
}

sub getSeqName {
  my $self = shift;
  return $self->{'seqname'};
}

1;
