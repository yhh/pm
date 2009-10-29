package entity;
use strict;
use Bio::Root::Root;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'entity';
$VERSION  = 0.1;

=head1 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($name, $type, $modifications) = 
      $self->_rearrange([qw(NAME TYPE MODIFICATIONS)], @args);
  
  $self->{'type'} = $type;
  $self->{'name'} = $name;
  $self->{'modifications'} = $modifications;

  return $self; # success - we hope!
}

sub _validate {
  #validate the startup parameters here.
}

sub getModifications {
  my $self  = shift;
  my @modifications = @{$self->{'modifications'}};
  return \@modifications;
}

sub getName {
  my $self  = shift;
  return $self->{'name'};
}

sub getType {
  my $self  = shift;
  return $self->{'type'};
}

1;
