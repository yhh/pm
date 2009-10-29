package simple_fasta;

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use genomedbrange;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root, Bio::Root::IO);
$ID = 'simple_fasta';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($assembly) =
     $self->_rearrange([qw(ASSEMBLY)], @args);

  $self->_initialize_io(@args);

  $self->{'assembly'} = $assembly;
  return $self; # success - we hope!
}

sub next_sequence {
   my ($self) = @_;
   while (my $aline=$self->_readline()) {
      chomp($aline); 
      if ($aline eq "") {next;}
      if ($aline=~/^>(\d+)/) {
         return $1;
      }
   }

   return undef;
}

sub getDBRList {
  my ($self) = @_;
  my @list=();

  while (my $aline=$self->_readline()) {
    chomp($aline);
    if ($aline eq "") {
       next;
    }

    if ($aline=~/^>/) {
       $self->_pushback($aline."\n");
       last;
    }

    my @fields = split("\t", $aline);
    my $aDBR = new genomedbrange(-DBVERSION => $self->{'assembly'},
                                 -SEQREGION => $fields[0],
                                 -POS1      => $fields[1],
                                 -POS2      => $fields[2],
                                 -STRAND    => $fields[3]
                                );

    if (!defined($aDBR)) {
       die "Something wrong with $aline\n";
    }

    push(@list, $aDBR);
  }
  return \@list;
}

1;

