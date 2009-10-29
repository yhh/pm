package tRNAClusterDB;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use tRNASeq;
use tRNASeqFasta;
use tRNAInfo;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root Bio::Root::IO);
$ID = 'tRNAClusterDB';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->_initialize_io(@args);
  my %clusterDB = ();
  my %memberDB  = ();
  my @order     = ();

  my $orderIdx = 0;
  while ($self->_hasNext()) {
     my $aSeq = $self->_nextSequence();
     for (my $i = 0; $i < $aSeq->length; $i++) {
        my $aSymbol = $aSeq->symbolAt($i + 1);
        $memberDB{$aSymbol->hashCode()} = $aSymbol;
     }
     
     $clusterDB{$aSeq->seqName} = $aSeq;
     push(@order, $aSeq->seqName);
  }
  
  $self->{'memberDB'}  = \%memberDB;
  $self->{'clusterDB'} = \%clusterDB;
  $self->{'order'}     = \@order;
  return $self; # success - we hope!
}

sub memberDB {
  my $self = shift;
  my $refmemberDB = $self->{'memberDB'};
  return $refmemberDB;
}

sub NCMemberDB {
  my $self = shift;
  my $refNCMemberDB = $self->{'NCMemberDB'};
  return $refNCMemberDB;
}

sub orderDB {
  my $self = shift;
  return $self->{'order'};
}

sub clusterDB {
  my $self = shift;
  my $refclusterDB = $self->{'clusterDB'};
  return $refclusterDB;
}

sub nonClusteredDB {
  my $self = shift;
  my $refNCDB = $self->{'nonClusteredDB'};
  return $refNCDB;
}

sub _hasNext {
  my $self = shift;
  while (my $aline = $self->_readline()) {
     chomp($aline);
     if ($aline =~/^>/) {
        $self->_pushback($aline."\n");
        return 1;
     }
  }
  return 0;
}

sub _nextSequence {
  my $self = shift;
  my $aSeq = undef;

  while (my $aline = $self->_readline()) {
     chomp($aline);
     if ($aline =~/^>(.*)/) {
        if (defined($aSeq)) {$self->_pushback($aline."\n"); last;}
        my $currentSeqName = join('.', split("\t", $1));
        $aSeq   = new tRNASeq(-symbollen=>5, -seqname=>$currentSeqName);
     } elsif ($aline eq '') {
        next;
     } elsif ($aline=~/^\!/) {
        my @fields = split("\t", $aline);
        shift @fields;

        my $ncName = shift @fields;

        #AA ANTICODON AALABEL CHR POS1 POS2 STRAND BITSCORE
        #Ser     CGA     Ser3F   12      54870415        54870496        1       89.14
        my $aSymbol = new tRNAInfo(-aa        => $fields[0],
                                   -anticodon => $fields[1],
                                   -aalabel   => $fields[2],
                                   '-chr'     => $fields[3],
                                   -pos1      => $fields[4],
                                   -pos2      => $fields[5],
                                   -strand    => $fields[6],
                                   -bitscore  => $fields[7]
                                  );


        $self->{'nonClusteredDB'}->{$ncName} = $aSymbol;
        $self->{'NCMemberDB'}->{$aSymbol->hashCode} = $aSymbol;

     } elsif ($aline=~/^\#\#(.*)/) {
        $aSeq->{'desc1'} = $1;
     } elsif ($aline=~/^\#(.*)/) {
	$aSeq->{'desc0'} = $1;
     } elsif ($aline!~/^\#/) {
        my @fields = split("\t", $aline);
        #pop @fields;
        $aSeq->addSymbol(\@fields);
     }
  }
  return $aSeq;
}

1;
