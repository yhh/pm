package eufindtRNA;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use tRNAHit;
use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root Bio::Root::IO);
$ID = 'eufindtRNA';
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

sub next_sequence {
  my $self = shift;
  while (my $aline = $self->_readline()) {
     chomp($aline);
     if ($aline=~/^Seq:\s+(.*?)$/) {
        $self->{'current_seqname'} = $1;
        return $self->{'current_seqname'};
     }
  }
  return undef;
}

sub next_hit {
  my $self = shift;
  my $seqname = $self->{'current_seqname'};
  $seqname=~s/\|/\\\|/g;
  my $aHit = undef;

  while (my $aline = $self->_readline()) {
    chomp($aline);
    if ($aline=~/^Bbox/) {
       if (!defined($aHit)) {
          $aHit = new tRNAHit(-seqname=>$self->{'current_seqname'});
       } else {
          $self->_pushback($aline."\n");
          last;
       }

       if (!$aHit->addBBox($aline)) {
           die "Something wrong with the output of eufindtRNA\n$aline\n";
       }
    } elsif ($aline=~/^Abox/) {
        #if (!defined($aHit)) {$aHit = new tRNAHit(-seqname=>$self->{'current_seqname'});}
        if (!$aHit->addABox($aline)) {
           die "Something wrong with the output of eufindtRNA\n$aline\n";
        }
    }

    if ($aline=~/^$seqname/) {
       $aHit->IsGoodHit(1);
       last;
    }
    
    if ($aline=~/^Seq:/) {
       $self->_pushback($aline."\n");
       if (defined($aHit)) {
          $aHit->IsGoodHit(0);
       }
       last;
    }
  }

  return $aHit;
}

1;
