package RawtRNASeqFasta;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use tRNASeq;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root Bio::Root::IO);
$ID = 'RawtRNASeqFasta';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->_initialize_io(@args);
  return $self; # success - we hope!
}

sub hasNext {
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

sub nextSequence {
  my $self   = shift;
  my $seqAln = {};
  my ($seqName1, $seqName2) = (undef, undef);
  my $seq1   = undef;
  my $seq2   = undef;

  while (my $aline = $self->_readline()) {
     chomp($aline);
     if ($aline eq '') {
        next;
     } elsif ($aline=~/^>(.*)/) {
        if (defined($seq1) && defined($seq2)) {
           $self->_pushback($aline."\n");
           last;
        }

        my @fields  = split("\t", $1);
        $seqName1 = join('.', $fields[0], $fields[1], $fields[2]);
        $seqName2 = join('.', $fields[4], $fields[5], $fields[6]);
        $seq1 = new tRNASeq(-symbollen=>5 , -seqname=>$seqName1);
        $seq2 = new tRNASeq(-symbollen=>5 , -seqname=>$seqName2);
        $seqAln->{'1'} = $seq1;
        $seqAln->{'2'} = $seq2;
     } else {
        my ($tempSeqName1, $pos1_1, $seqString1, $pos1_2) = split(/\s+/, $aline);
        
        my $aline_2 = $self->_readline();
        chomp($aline_2);
        my ($tempSeqName2, $pos2_1, $seqString2, $pos2_2) = split(/\s+/, $aline_2);

        $self->addSymbols($seqString1, $seq1);
        $self->addSymbols($seqString2, $seq2);
     }
  }
  return $seqAln;
}

sub addSymbols {
  my $self = shift;
  my $seqString  = shift;
  my $tRNASeqObj = shift;
  while ($seqString =~/(.{4})(.)/g) {
     my @params = ($1, undef, $1.$2);
     $tRNASeqObj->addSymbol(\@params);
  }
}

1;
