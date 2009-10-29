package tRNASeq;
use strict;
use Bio::Root::Root;
use tRNAInfo;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'tRNASeq';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  ($self->{'symbolLen'}, $self->{'seqName'})= $self->_rearrange([qw(SYMBOLLEN SEQNAME)], @args);
  if (!defined($self->{'symbolLen'})) {$self->{'symbolLen'} = 5;}
  if (!defined($self->{'seqName'})) {return undef;}
  
  ($self->{'clusterID'}, $self->{'chr'}, $self->{'memberCount'}) = split('\.', $self->{'seqName'});
  $self->{'seqString'} = '';
  return $self; # success - we hope!
}

sub seqString {
  my $self = shift;
  my $seqString = shift;
  if (defined($seqString)) {
     $self->{'seqString'} = $seqString;
  } else {
     return $self->{'seqString'};
  }
}

sub seqName {
  my ($self, $name) = @_;
  if (defined($name)) {
     $self->{'seqName'} = $name;
  } else {
     return $self->{'seqName'};
  }
}

sub addSymbol {
  my ($self, $refInf) = @_;
  #print $refInf->[2], "\n";
  my $aSymbol = new tRNAInfo(-aa        => $refInf->[0],     
                             -anticodon => $refInf->[1], 
                             -aalabel   => $refInf->[2], 
                             "-chr"     => $refInf->[3],  
                             -pos1      => $refInf->[4],
                             -pos2      => $refInf->[5], 
                             -strand    => $refInf->[6], 
                             -bitscore  => $refInf->[7]);
  if ($self->{'symbolLen'} != length($refInf->[2])) {
     warn "Symbol length is not compatible with ".$refInf->[2]." \n";
     return;
  }

  $self->{'seqString'}.= $refInf->[2];
  #print $aSymbol, "\n";
  push(@{$self->{'seqSyms'}}, $aSymbol);
  $self->{'length'}++;
  if ($aSymbol->{'aa'} !~/\-{2,}/) {
     #print $aSymbol->{'aa'}, "\n";
     $self->{'real_length'}++;
  } 
}

sub firstSymbol {
  my $self = shift;
  return $self->symbolAt(1);
}

sub lastSymbol {
  my $self = shift;
  return $self->symbolAt($self->length);
}

sub dist {
  my $self  = shift;
  my $idx   = shift;
  my $isRev = shift;

  if ($isRev) {
     my $next_idx = $idx + 1;
     if ($next_idx > $self->length) {
        return 0;
     }
     
     return $self->symbolAt($next_idx)->pos1 - $self->symbolAt($idx)->pos2;
  } else {
     my $prev_idx = $idx - 1;
     if ($prev_idx < 1) {
        return 0;
     }

     return $self->symbolAt($idx)->pos1 - $self->symbolAt($prev_idx)->pos2;
  }
}

sub symbolAt {
  my ($self, $idx, $isReverse) = @_;
  if ($isReverse) {
     $idx = $self->length - $idx +1;
  }

  if ($idx > $self->length || $idx < 1) {
     return undef;
  }
  
  return $self->{'seqSyms'}->[$idx-1];
}

sub withinSeqRange {
  my ($self, $client_chr, $client_coord1) = shift;
  my $firstSymbol = $self->firstSymbol;
  my $lastSymbol  = $self->lastSymbol;
  my $chr         = $self->firstSymbol->chr;
  
  if ($client_chr ne $chr) {
     return 0;
  }

  my @coords = sort ($firstSymbol->pos1, $firstSymbol->pos2, 
                     $lastSymbol->pos1,  $lastSymbol->pos2);
  
  if ($client_coord1 >= $coords[0] && $client_coord1 <= $coords[3]) {
     return 1;
  } else {
     return 0;
  }
}

sub DBRange {
  my $self = shift;
  my $assembly = shift;
  
  my $firstDBR = $self->firstSymbol->GDBR($assembly);
  my $lastDBR  = $self->lastSymbol->GDBR($assembly);

  my $mergedDBR = $firstDBR->mergeDBR($lastDBR);
  return $mergedDBR;
}

sub labelList {
  my $self = shift;
  my @arrayOfLabel = ();
  for (my $i=0; $i < scalar(@{$self->{'seqSyms'}}); $i++) {
      push(@arrayOfLabel, $self->{'seqSyms'}->[$i]->aalabel);
  }
  return \@arrayOfLabel;
}

sub length {
  my $self = shift;
  return $self->{'length'};
}

sub real_length {
  my $self = shift;
  return $self->{'real_length'};
}

sub getChr {
  my $self = shift;
  return $self->{'chr'};
}

sub getClusterID {
  my $self = shift;
  return $self->{'clusterID'};
}

sub getClusterMemberCount {
  my $self = shift;
  return $self->{'memberCount'};
}

1;
