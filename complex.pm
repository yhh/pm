package complex;
use strict;
use Bio::Root::Root;
use entity;
use Clone::Fast qw( clone );

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'complex';
$VERSION  = 0.1;

=head1 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($complex_ID, $entities) = 
      $self->_rearrange([qw(COMPLEX_ID ENTITIES)], @args);
  
  $self->{'complex_id'} = $complex_ID;
  $self->{'entities'}   = $entities;

  return $self; # success - we hope!
}

sub getEntities {
  my $self = shift;
  my $entities = $self->{'entities'};
  my $clones   = [];
  foreach my $single_entity (@$entities) {
     push(@$clones, clone($single_entity));
  }
  return $clones;
}

sub _validate {
  #validate startup parameters here.
}

sub getComplexID {
  my $self = shift;
  return $self->{'complex_id'};
}

1;
