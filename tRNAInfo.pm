package tRNAInfo;

use strict;
use Bio::Root::Root;
use genomedbrange;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'tRNAInfo';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  ($self->{'aa'}, $self->{'anticodon'}, $self->{'aalabel'}, 
   $self->{'chr'},$self->{'pos1'},      $self->{'pos2'}, $self->{'strand'}, $self->{'bitscore'}) = 
     $self->_rearrange([qw(AA ANTICODON AALABEL CHR POS1 POS2 STRAND BITSCORE)], @args);
  return $self; # success - we hope!
}

sub aa {
  my $self = shift;
  return $self->{'aa'};
}

sub hashCode {
  my $self = shift;
  if (!defined($self->{'hashcode'})) {
     $self->{'hashcode'}=join('.', $self->chr, $self->pos1, $self->pos2, $self->strand);
  }
  return $self->{'hashcode'};
}

sub GDBR {
  my $self = shift;
  my $assembly = shift;
  if ($assembly eq "") {
     warn "no genome assembly provided!\n";
     return undef;
  } else {
     #SEQLEVEL DBVERSION SEQREGION POS1 POS2 STRAND)
     return new genomedbrange(-SEQLEVEL  => 'chromosome', 
                              -DBVERSION => $assembly, 
                              -SEQREGION => $self->chr, 
                              -POS1      => $self->pos1,
                              -POS2      => $self->pos2,
                              -STRAND    => $self->strand);
  }
}

sub anticodon {
  my $self = shift;
  return $self->{'anticodon'};
}

sub aalabel {
  my $self = shift;
  return $self->{'aalabel'};
}

sub aalabel_trim {
  my $self = shift;
  my $aalabel_trim = $self->{'aalabel'};
  chop($aalabel_trim);
  return $aalabel_trim;
}

sub chr {
  my $self = shift;
  return $self->{'chr'};
}

sub pos1 {
  my $self = shift;
  return $self->{'pos1'};
}

sub pos2 {
  my $self = shift;
  return $self->{'pos2'};
}

sub strand {
  my $self = shift;
  return $self->{'strand'};
}

sub bitscore {
  my $self = shift;
  return $self->{'bitscore'};
}

sub toString {
  my $self = shift;
  return join("\t", $self->chr, $self->pos1, $self->pos2, $self->strand, $self->bitscore, $self->aalabel, $self->anticodon);
}

sub toString2 {
  my $self = shift;
  return join("\t", $self->aa, 
                    $self->anticodon,
                    $self->aalabel,
                    $self->chr, 
                    $self->pos1, 
                    $self->pos2,
                    $self->strand,
                    $self->bitscore);
}

1;
