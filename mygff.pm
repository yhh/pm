package mygff;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use simplerange;
use vars qw($ID $VERSION @ISA $filter);

@ISA= qw(Bio::Root::Root Bio::Root::IO);
$ID = 'mygff';
$filter = undef;
$VERSION  = 0.1;

=head1 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # initialize IO
  $self->_initialize_io(@args);

  return $self; # success - we hope!
}

sub hasNext {
   my ($self) = @_;
   while (my $aline=$self->_readline()) {
     chomp($aline);
     if ($aline ne "") {
        if (defined($filter) && $aline!~/$filter/) {
           next;
        }

        if ($aline=~/^#/) {
           next;
        } else {
           $self->_pushback($aline."\n");
           return 1;
        }
     }
  }
  return 0;
}

sub setFilter {
  my $self;
  ($self, $filter) = @_;
}

sub nextRange {
  my ($self) = @_;

  my $aline = $self->_readline();
  chomp($aline);

  my ($name, $program, $model, $start, $end, $prob, $strand, $desc) = split("\t", $aline);
  if ($name=~/\/(\d+)-(\d+)/) {
     $start+=$1-1;
     $end+=$1-1;
  }

  if ($strand eq '-') {
     $strand = -1;
  } else {
     $strand = 1;
  }

  my $aRange = new simplerange(-pos1=>$start, -pos2=>$end, -strand=>$strand, -geneclass=>$name, -desc=>$aline);

  return $aRange;
}

1;
