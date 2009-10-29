package alignedSymbols;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root Bio::Root::IO);
$ID = 'tRNAClusterDB';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);  
  $self->_initialize_io(@args);
  return $self;
}

sub _hasNext {
  my $self = shift;
  while (my $aline = $self->_readline()) {
     chomp($aline);
     if ($aline=~/^>/) {
        $self->_pushback($aline."\n");
        return 1;
     }
  }
  return 0;
}

sub _nextAlignment {
  my $self = shift;
  my $aAln = undef;

  while (my $aline = $self->_readline()) {
     chomp($aline);
     if ($aline =~/^>(.*?)\s+<=>\s+(.*)/) {
		if (defined($aAln)) {$self->_pushback($aline."\n"); last;}

		my $name1 = $1;
		my $name2 = $2;

		$name1=~s/\s+/\./g;
		$name2=~s/\s+/\./g;

		$aAln = {};

		if ($name2=~/(.*?)\_R/) {
			$name2 = $1;
			$aAln->{'isRev'} = 1;
		} else {
			$aAln->{'isRev'} = 0;
		}
		$aAln->{'name1'} = $name1;
		$aAln->{'name2'} = $name2;

		$aAln->{'seq1'} = ();
		$aAln->{'seq2'} = ();
		$aAln->{'seq1_order'} = ();
		$aAln->{'seq2_order'} = ();

	 } elsif ($aline ne "") {
		my $aline1  = $self->_readline();
		my $aline2 = $self->_readline();
		my $dist2 = $self->_readline();
		
		chomp($aline1);
		chomp($aline2);
		chomp($dist2);
		
		#print $aline, "\n";
		#print $aline1, "\n";

		if ($aline1=~/\S+\s+\d+\s+(.*?)\s+\d+/) {
			push(@{$aAln->{'seq1'}}, split('~', $1));
		}

		if ($aline2=~/\S+\s+\d+\s+(.*?)\s+\d+/) {
			push(@{$aAln->{'seq2'}}, split('~', $1));
		}
	 }
  }

  my $seq1_idx = 0;
  my @seq1_order = ();
  foreach my $single_symbol (@{$aAln->{'seq1'}}) {
     if ($single_symbol !~/\-\-\-/) {
        $seq1_idx++;
        push(@seq1_order, $seq1_idx);
     } else {
        push(@seq1_order, -1);
     }
  }
  $aAln->{'seq1_order'} = \@seq1_order;
  $aAln->{'seq1_len'}   = $seq1_idx;

  my $seq2_idx = 0;
  my @seq2_order = ();
  foreach my $single_symbol (@{$aAln->{'seq2'}}) {
     if ($single_symbol !~/\-\-\-/) {
        $seq2_idx++;
        push(@seq2_order, $seq2_idx);
     } else {
        push(@seq2_order, -1);
     }
  }
  $aAln->{'seq2_len'}   = $seq2_idx;

  if ($aAln->{'isRev'}) {
     my @rev_seq2_order = ();
     foreach my $single_value (@seq2_order) {
        if ($single_value < 0) {
           push(@rev_seq2_order, $single_value);
        } else {
           push(@rev_seq2_order, $seq2_idx - $single_value + 1);
        }
     }
     $aAln->{'seq2_order'} = \@rev_seq2_order;
  } else {
     $aAln->{'seq2_order'} = \@seq2_order;
  }

  return $aAln;
}

1;
