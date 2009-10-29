package tRNASymbol;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'infernal';
$VERSION  = 0.1;

my $sequence_name=undef;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # initialize IO
  $self->_initialize_io(@args);

  return $self; # success - we hope!
}




1;

