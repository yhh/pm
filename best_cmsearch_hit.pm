package best_cmsearch_hit;

use strict;
use genomedbrange;
use Bio::Root::Root;
use infernal3;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'best_cmsearch_hit';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my $if   = new infernal3(@args);

  my @hits = ();
  my $hit_selected;

  while ($if->next_sequence) {
    while (my $hit=$if->next_hit) {
      my $dbrange = new genomedbrange(-dbrange => $hit->{'sequence_name'});

      my $stranded_start  = $hit->{'hit_start'};
      my $stranded_end    = $hit->{'hit_end'};
      my $stranded_strand = 1;
      if ($stranded_start > $stranded_end) {
         $stranded_strand = -1;
         my $temp = $stranded_start;
         $stranded_start = $stranded_end;
         $stranded_end   = $temp;
      }

      #$dbrange    = $dbrange->getSubRange($hit->{'hit_start'}, $hit->{'hit_end'});
      #$dbrange    = $dbrange->getSubRange($stranded_start, $stranded_end, $stranded_strand);
      #$hit->{'genomedbrange'} = $dbrange;
      push(@hits, $hit);
      if (!defined($hit_selected)) {
         $hit_selected = $hit;
      } else {
         $hit_selected = ($hit_selected->{'bits'} > $hit->{'bits'} 
                        ? $hit_selected           : $hit);
      }
    }
  }

  @hits = sort {$b->{'bits'} <=> $a->{'bits'}} @hits;

  $self->{'hitcount'}         = scalar(@hits);
  $self->{'hit_selected'}     = $hit_selected;
  $self->{'hits'}             = \@hits;
  return $self;
}

sub toString {
  my $self = shift;
  my $hitName;
  if (!defined($self->getBestHitGenomeDBRange)) {
     $hitName = "undefRfam";
  } else {
     $hitName = $self->getBestHitGenomeDBRange->getSequenceName;
  }
  return join("\t", $self->numHits,
                    $self->bestBits,
                    $hitName
             );
}

sub hitSelected {
  my $self = shift;
  return $self->{'hit_selected'};
}

sub numHits {
  my $self = shift;
  if (!defined($self->{'hitcount'})) {
     return 0;
  } else {
     return $self->{'hitcount'};
  }
}

sub bestBits {
  my $self = shift;
  if (!defined($self->hitSelected)) {
     return 0;
  } else {
     return $self->hitSelected->{'bits'};
  }
}

sub getBestHitGenomeDBRange {
  my $self = shift;
  if (!defined($self->hitSelected)) {
     return undef;
  } else {
     return $self->hitSelected->{'genomedbrange'};
  }
}

1;
