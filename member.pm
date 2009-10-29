package member;
use strict;
use Bio::Root::Root;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'member';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($id, $score, $desc) = $self->_rearrange([qw(ID SCORE DESC)], @args);
  $self->{'id'}    = $id;
  $self->{'score'} = $score;
  $self->{'desc'}  = $desc;
  return $self;
}

sub id {
  my $self = shift;
  return $self->{'id'};
}

sub score {
  my $self = shift;
  return $self->{'score'};
}

sub desc {
  my $self = shift;
  return $self->{'desc'};
}

1;
