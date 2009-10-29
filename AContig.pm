package AContig;

use strict;
use Bio::Root::Root;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'AContig';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($contig) = 
     $self->_rearrange([qw(CONTIG)], @args);
  
  if ($contig=~/(.*?)\.(\d+)\.(\d+)\.(\d+)/) {
     $self->{'contig_acc'} = $1.'.'.$2;
     $self->{'start'}  = $3;
     $self->{'end'}    = $4;
     $self->{'strand'} = 1;

     if ($self->{'start'} > $self->{'end'}) {
        my $temp = $self->{'start'};

        $self->{'start'} = $self->{'end'};
        $self->{'end'}   = $temp;

        $self->{'strand'} = -1;
     }
  } else {
     warn "$contig is not a valid range\n";
     $self = undef;
  }

  return $self; # success - we hope!
}

sub new2 {
  my ($class, @args) =@_;
  my $self = $class->SUPER::new(@args);
 
  my ($contig_acc, $start, $end, $strand) = 
     $self->_rearrange([qw(CONTIG_ACC START END STRAND)], @args); 
  
  if (defined($contig_acc) && defined($start) && defined($end)) {

     $self->{'contig_acc'} = $contig_acc;
     $self->{'start'}      = $start;
     $self->{'end'}        = $end;
     $self->{'strand'}     = 1;
     
     
     if ($self->{'start'} > $self->{'end'}) {
        my $temp = $self->{'start'};

        $self->{'start'} = $self->{'end'};
        $self->{'end'}   = $temp;

        $self->{'strand'} = -1;

     } elsif (defined($strand)) {
        $self->{'strand'} = $strand;
     }
  } else {
     return undef;
  }
  
  return $self;
}

sub contig_acc {
  my $self = shift;
  return $self->{'contig_acc'};
}

sub start {
  my $self = shift;
  return $self->{'start'};
}

sub end {
  my $self = shift;
  return $self->{'end'};
}

sub strand {
  my $self = shift;
  return $self->{'strand'};
}

sub toString {
  my $self = shift;
  return $self->contig_acc.'.'.$self->start.'.'.$self->end.'.'.$self->strand;
}

sub subRange {
  my $self = shift;
  my ($pos1, $pos2) = (shift, shift);
  
  if (!defined($pos1) || !defined($pos2)) {
     die "subrange pos1 or pos2 is wrong in AContig!\n";
  }

  my $new_start = $pos1 + $self->start - 1;
  my $new_end   = $pos2 + $self->start - 1;

  if ($new_end > $self->end) {
     warn "$new_end is behind old end".$self->end."\n";
     $new_end = $self->end;
  }

  return $self->new2(-contig_acc=>$self->contig_acc, -start=>$new_start, -end=>$new_end);
}
