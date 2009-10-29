package tRNASeqFasta;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use tRNASeq;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root Bio::Root::IO);
$ID = 'tRNASeqFasta';
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
        my ($seqName1, $seqName2)  = split("\t", $1);
        #$seqName1 = join('.', $fields[0], $fields[1], $fields[2]);
        #$seqName2 = join('.', $fields[4], $fields[5], $fields[6]);
        $seq1 = new tRNASeq(-symbollen=>5 , -seqname=>$seqName1);
        $seq2 = new tRNASeq(-symbollen=>5 , -seqname=>$seqName2);
        $seqAln->{'1'} = $seq1;
        $seqAln->{'2'} = $seq2;
     } else {
        #my ($tempSeqName1, $pos1_1, $seqString1, $pos1_2) = split(/\s+/, $aline);
        my $seqString1 = $aline;
        my $aline_2 = $self->_readline();
        chomp($aline_2);
        #my ($tempSeqName2, $pos2_1, $seqString2, $pos2_2) = split(/\s+/, $aline);
        my $seqString2 = $aline_2;
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

sub identifyContinuousBlock {
  my $self = shift;
  my $aAln = shift;
  my $lengthcriteria = shift;
 
  if (!defined($lengthcriteria)) {
     $lengthcriteria = 2;
  }
 
  my @ranges = ();
  
  my @labelList1 = @{$aAln->{'1'}->labelList};
  my @labelList2 = @{$aAln->{'2'}->labelList};
  
  #print join("\t", @labelList1), "\n";
  #print join("\t", @labelList2), "\n";

  my $startIdx1;
  my $endIdx1;
  my $startIdx2;
  my $endIdx2;
  
  my $realIdx1;
  my $realIdx2;

  #print scalar(@labelList1), "\n";
  #print scalar(@labelList2), "\n";
  for (my $i = 0; $i < scalar(@labelList1); $i++) {
     if ($labelList1[$i]!~/\-\-\-/) {$realIdx1++;}
     if ($labelList2[$i]!~/\-\-\-/) {$realIdx2++;}
     
     if ($labelList1[$i]!~/\-\-\-/ && $labelList2[$i]!~/\-\-\-/) {
        #print join("\t", $i, $labelList1[$i], $labelList2[$i]), "\n";
        if (!defined($startIdx1)) {
           $startIdx1 = $realIdx1;
           $startIdx2 = $realIdx2;
        } else {
           $endIdx1   = $realIdx1;
           $endIdx2   = $realIdx2;
        }
     } elsif (defined($startIdx1)) {
        if ($endIdx1 - $startIdx1 + 1 >= $lengthcriteria) {
           my @temp = ($startIdx1, $endIdx1, $startIdx2, $endIdx2);
           push(@ranges, \@temp);
        }
        $startIdx1 = undef;
        $endIdx1   = undef;
        $startIdx2 = undef;
        $endIdx2   = undef;
     }
  }

  if (defined($endIdx1)) {
     if ($endIdx1 - $startIdx1 + 1 > $lengthcriteria) {
        my @temp = ($startIdx1, $endIdx1, $startIdx2, $endIdx2);
        push(@ranges, \@temp);
     }
  }
  return \@ranges;
}

1;
