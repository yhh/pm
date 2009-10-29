package checkoverlap;
use strict;
use Bio::Root::Root;
use math;
use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'checkoverlap';
$VERSION  = 0.1;

=head1 

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my $refRanges = $self->_rearrange( ["refranges"], @args);
  $self->{"refranges"} = $refRanges;
  return $self; # success - we hope!
}

sub isOverlap {
  my ($self, @args) = @_;
  my $aRange = $args[0];
  if (!$aRange->isa('genomedbrange')) {
     warn "Please check the range you've given is a genomedbrange\n";
     return 0;
  }
 
  my $overlaplength = 0;
  foreach my $aFeature (@{$self->getRefRanges}) {
     if (!$aFeature->isa("Bio::EnsEMBL::Feature")) {
        warn "Please check your features entered!\n";
        next;
     }
     
     #print join("\t", $aFeature->start, $aFeature->end), "\n";
     my $overlapstart = max($aFeature->start + $aRange->getPos1 -1, $aRange->getPos1);
     my $overlapend   = min($aFeature->end   + $aRange->getPos1 -1, $aRange->getPos2);
     
     $overlaplength = $overlapend - $overlapstart +1;
     
     if ($overlaplength >0) {
        return 1;
     }
  }
  return 0;
}

sub getRefRanges {
  my $self = shift;
  return $self->{'refranges'};
}

1;
