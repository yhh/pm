package best_paired_cmsearch_hit;

use strict;
use genomedbrange;
use Bio::Root::Root;
use best_cmsearch_hit;
use math;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'best_paired_cmsearch_hit';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($fileName1, $fileName2) = $self->_rearrange([qw(file1 file2)], @args);
  my $bch1      = new best_cmsearch_hit(-file=>$fileName1);
  my $bch2      = new best_cmsearch_hit(-file=>$fileName2);

  #blank result here
  #$self->{'info_count'} = 0;
  #$self->{'cm_count'}   = 0;
  #$self->{'pcm_count'}  = 0;
  #$self->{'mis_count'}  = 0;
  #using $self->{'struc_info'}->... instead.
  
  $self->{'cm1'} = $bch1;
  $self->{'cm2'} = $bch2;

  if ($bch1->numHits == 0 || $bch2->numHits == 0) {
     #blank result here
  } else {
     #trim insertions relative to the standard CM strucs.....
     my $res = doStrucAln($bch1, $bch2);
     $self->addTrimAln($res); #no ref will be returned, new refs are added directly to the original ref

     #put $res to somewhere in $self?
     
     #get CMs, update info_count, cm_count, pcm_count, and mis_count here, maybe conserved_insert_count here
     $self->getCMs($res);
  }

  return $self;
}

sub newByPairedObj {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($bch1, $bch2) = $self->_rearrange([qw(bch1 bch2)], @args);
  $self->{'cm1'} = $bch1;
  $self->{'cm2'} = $bch2;
  
  if ($bch1->numHits == 0 || $bch2->numHits == 0) {
     #blank result here
  } else {
     my $res = doStrucAln($bch1, $bch2);
     $self->addTrimAln($res);
     $self->getCMs($res);
  }
  return $self;
}

sub cm1Bits {
  my $self = shift;
  return $self->{'cm1'}->bestBits;
}

sub cm2Bits {
  my $self = shift;
  return $self->{'cm2'}->bestBits;
}

sub cmCount {
  my $self = shift;
  my $temp = $self->{'struc_info'}->{'CM'};
  if (!defined($temp)) {
     $temp = 0;
  }
  return $temp;
}

sub pcmCount {
  my $self = shift;
  my $temp =  $self->{'struc_info'}->{'partial_CM'};
  if (!defined($temp)) {
     $temp = 0;
  }
  return $temp;
}

sub stNoChange {
  my $self = shift;
  my $temp = $self->{'struc_info'}->{'st_no_change'};
  if (!defined($temp)) {
     $temp = 0; 
  }
  return $temp;
}

sub stMismatch {
  my $self = shift;
  my $temp = $self->{'struc_info'}->{'st_mismatch'};
  if (!defined($temp)) { 
     $temp = 0;
  }
  return $temp;
}

sub stInfoCount {
  my $self = shift;
  my $temp1 = $self->{'stem_len'};
  if (!defined($temp1)) {
     $temp1 = 0;
  }

  my $temp2 = $self->{'struc_info'}->{'st_no_info'};
  if (!defined($temp2)) {
     $temp2 = 0;
  }

  return $temp1 - $temp2;
}

sub stNoInfo {
  my $self = shift;
  my $temp = $self->{'struc_info'}->{'st_no_info'};
  if (!defined($temp)) {
     $temp = 0;
  }
  return $temp;
}

sub ssNoChange {
  my $self = shift;
  my $temp = $self->{'struc_info'}->{'ss_no_change'};
  if (!defined($temp)) {
     $temp = 0;
  }
  return $temp;
}

sub ssMismatch {
  my $self = shift;
  my $temp = $self->{'struc_info'}->{'ss_mismatch'};
  if (!defined($temp)) {
     $temp = 0;
  }
  return $temp;
}

sub ssInfoCount {
  my $self = shift;

  my $temp1 = $self->{'ss_len'};
  my $temp2 = $self->{'struc_info'}->{'ss_no_info'};
  if (!defined($temp1)) {
     $temp1 = 0;
  }

  if (!defined($temp2)) {
     $temp2 = 0;
  }
  return $temp1 - $temp2;
}

sub ssNoInfo {
  my $self = shift;
  my $temp = $self->{'struc_info'}->{'ss_no_info'};
  if (!defined($temp)) {
     $temp = 0;
  }
  return $temp;
}

sub stLen {
  my $self = shift;
  if (!defined($self->{'stem_len'})) {
     return 0;
  }
  return $self->{'stem_len'};
}

sub ssLen {
  my $self = shift;
  if (!defined($self->{'ss_len'})) {
     return 0;
  }
  return $self->{'ss_len'};
}

sub alnLen {
  my $self = shift;
  return $self->{'aln_len'}
}

sub getCMs {
  my $self = shift;
  my $res  = shift;
  my %struc_info=();

  $self->{'struc_info'} = \%struc_info;

  #  {'st_no_info'}
  #  {'st_no_change'}
  #  {'st_mismatch'}
  #  {'partial_CM'}
  #  {'CM'}
  #  {'st_strange'}
  #  {'st_very_strange'}
  #  {'ss_no_info'}
  #  {'ss_no_change'}
  #  {'ss_mismatch'}
  #  checksum?=================> 
  #  2 * (st_no_info + st_no_change + st_mismatch + partial_CM + CM + st_strange + st_very_strange) + 
  #  ss_no_info + ss_no_change + ss_mismatch

  #get paired indices
  my %pairhashforward=();
  my %pairhashreverse=();
  for (my $j=0; $j<2; $j++) {
     my @struc_stack=();
     for (my $seqidx=0; $seqidx<scalar(@{$res->[$j]->{'tr_struc'}}); $seqidx++) {
         my $temp_struc=$res->[$j]->{'tr_struc'}->[$seqidx];
         if ($temp_struc eq '<' || $temp_struc eq '(') {
            push(@struc_stack, $seqidx);
         } elsif ($temp_struc eq '>' || $temp_struc eq ')') {
            my $seqidx0 = pop @struc_stack;
            #my $paired_coords=$seqidx0."\t".$seqidx;
            #print $paired_coords, "\n";
            $pairhashforward{$seqidx0}=$seqidx;
            $pairhashreverse{$seqidx} =$seqidx0;
         }
     }
  }

  $self->{'stem_len'} = scalar(keys %pairhashforward);
  $self->{'aln_len'}  = scalar(@{$res->[0]->{'tr_struc'}});
  $self->{'ss_len'}   = $self->{'aln_len'} - $self->{'stem_len'} * 2;

  for (my $seqidx = 0; $seqidx < scalar(@{$res->[0]->{'tr_struc'}}); $seqidx++) {
     my $pair1_1 = $res->[0]->{'tr_targetseq'}->[$seqidx];
     my $pair2_1 = $res->[1]->{'tr_targetseq'}->[$seqidx];

     if (defined($pairhashforward{$seqidx}) || defined($pairhashreverse{$seqidx})) {
        #get information about the double stranded regions
        if (defined($pairhashforward{$seqidx})) {
           #get CMs
           #get paired bases
           #my $pair1_1 = $res->[0]->{'tr_targetseq'}->[$seqidx];
           my $pair1_2 = $res->[0]->{'tr_targetseq'}->[$pairhashforward{$seqidx}];
           #my $pair2_1 = $res->[1]->{'tr_targetseq'}->[$seqidx];
           my $pair2_2 = $res->[1]->{'tr_targetseq'}->[$pairhashforward{$seqidx}];

           #check no_info
           if ($pair1_1=~/[\s\-]/ || $pair1_2=~/[\s\-]/  || $pair2_1=~/[\s\-]/ || $pair2_2=~/[\s\-]/) {
              $struc_info{'st_no_info'}++;
           } elsif ($pair1_1 eq $pair2_1 && $pair1_2 eq $pair2_2) {
              $struc_info{'st_no_change'}++;
           } elsif ($pair1_1 ne $pair2_1 || $pair1_2 ne $pair2_2) {
              if ((!IsPaired($pair1_1, $pair1_2)) || (!IsPaired($pair2_1, $pair2_2))) {
                 $struc_info{'st_mismatch'}++;
              } elsif ($pair1_1 eq $pair2_1 || $pair1_2 eq $pair2_2) {
                 $struc_info{'partial_CM'}++;
              } elsif ($pair1_1 ne $pair2_1 && $pair1_2 ne $pair2_2) {
                 $struc_info{'CM'}++;
              } else {
                 $struc_info{'st_strange'}++;
              }
           } else {
              $struc_info{'st_very strange'}++;
           }
        } # in stem regions
     } else {
        #get information about the single stranded regions
        if ($pair1_1=~/[\s\-]/ || $pair2_1=~/[\s\-]/) {
           $struc_info{'ss_no_info'}++;
        } elsif ($pair1_1 eq $pair2_1) {
           $struc_info{'ss_no_change'}++;
        } else {
           $struc_info{'ss_mismatch'}++;
        }
     } # in single stranded regions
  } # iterate over all bases

  # process inserts or deletions?
  

  return \%struc_info;
}

sub addTrimAln {
  my $self   = shift;
  my $orires = shift;

     #$res->[$i]->{'stdseqtexts'} = $stdseqtexts[$i];
     #$res->[$i]->{'structexts'}  = $structexts[$i];
     #$res->[$i]->{'CMtexts'}     = $CMtexts[$i];
     #$res->[$i]->{'seqtexts'}    = $seqtexts[$i];
  
  my $oriseq= [{}, {}];
  my $trimed_seqs = [{}, {}];

  for (my $i = 0; $i < 2; $i++) {
      my @temp =split(//, $orires->[$i]->{'structexts'});   $oriseq->[$i]->{'struc'}    = \@temp;
      my @temp2=split(//, $orires->[$i]->{'stdseqtexts'});  $oriseq->[$i]->{'stdseq'}   = \@temp2;
      my @temp3=split(//, $orires->[$i]->{'CMtexts'});      $oriseq->[$i]->{'CMstring'} = \@temp3;
      my @temp4=split(//, $orires->[$i]->{'seqtexts'});     $oriseq->[$i]->{'targetseq'}= \@temp4;
      #print @temp, "\n";

      $trimed_seqs->[$i] = trim_extra_seq($oriseq->[$i]);

      $orires->[$i]->{'tr_struc'}       = $trimed_seqs->[$i]->{'struc'};
      $orires->[$i]->{'tr_stdseq'}      = $trimed_seqs->[$i]->{'stdseq'};
      $orires->[$i]->{'tr_CMstring'}    = $trimed_seqs->[$i]->{'CMstring'};
      $orires->[$i]->{'tr_targetseq'}   = $trimed_seqs->[$i]->{'targetseq'};
      $orires->[$i]->{'tr_inserts_ref'} = $trimed_seqs->[$i]->{'inserts_ref'};
  }
  
  #validate struc here.....
  my $len_for_checking = min(scalar(@{$trimed_seqs->[0]->{'struc'}}), 
                             scalar(@{$trimed_seqs->[1]->{'struc'}})
                            );
  my $validStrucAln = "validStrucAln";  

  for (my $i = 0; $i < $len_for_checking; $i++) {
      my $struc0 = $trimed_seqs->[0]->{'struc'}->[$i];
      my $struc1 = $trimed_seqs->[1]->{'struc'}->[$i];
      if ($struc0 ne " " && $struc1 ne " ") {
         if ($struc0 ne $struc1) {
            $validStrucAln = "warnedStrucAln";
            last;
         }
      }
  }

  $self->{'validStrucAln'} = $validStrucAln;

  return $orires;
}

sub trim_extra_seq {
  my $seq=shift;
  my $newseq;
  my $pseudo_idx=-1;
  my $inserts_ref = {};

  for (my $i=0; $i<scalar(@{$seq->{'struc'}}); $i++) {
      if ($seq->{'struc'}->[$i] ne '.') {
         $pseudo_idx++;
         $newseq->{'struc'}->[$pseudo_idx]     = $seq->{'struc'}->[$i];
         $newseq->{'stdseq'}->[$pseudo_idx]    = $seq->{'stdseq'}->[$i];
         $newseq->{'CMstring'}->[$pseudo_idx]  = $seq->{'CMstring'}->[$i];
         $newseq->{'targetseq'}->[$pseudo_idx] = $seq->{'targetseq'}->[$i];
      } else {
         $inserts_ref->{$i} = $seq->{'targetseq'}->[$i];
      }
  }
  $newseq->{'inserts_ref'} = $inserts_ref;
  return $newseq;
}

sub IsPaired {
  my $base1=uc(shift);
  my $base2=uc(shift);
  my $RC_base2=$base2;
  $RC_base2=~tr/AUTGC/UAACG/;
  if ($base1 eq $RC_base2) {
     return 1;
  } elsif ($base1 eq 'G') {
     if ($RC_base2 eq 'A') {
        return 1;
     } else {
        return 0;
     }
  } elsif ($base1 eq 'U') {
     if ($RC_base2 eq 'C') {
       return 1;
     } else {
        return 0;
     }
  } else {
     return 0;
  }
}


sub doStrucAln {
  my $bch1 = shift;
  my $bch2 = shift;
  my $res  = [{}, {}];

  my @best_hits = ($bch1->hitSelected, $bch2->hitSelected);
  
  my $supposed_seq_end  =($best_hits[0]->{'std_end'}   >  $best_hits[1]->{'std_end'} 
                        ? $best_hits[0]->{'std_end'}   :  $best_hits[1]->{'std_end'});

  my $supposed_seq_start=($best_hits[0]->{'std_start'} <  $best_hits[1]->{'std_start'}
                        ? $best_hits[0]->{'std_start'} :  $best_hits[1]->{'std_start'});

  my @oriseqs   =($best_hits[0]->{'target_seq'}, $best_hits[1]->{'target_seq'});
  my @oristdseqs=($best_hits[0]->{'std_seq'},    $best_hits[1]->{'std_seq'});
  my @oristrucs =($best_hits[0]->{'struc'},      $best_hits[1]->{'struc'});
  my @oriCMs    =($best_hits[0]->{'CMstring'},   $best_hits[1]->{'CMstring'});
  
  my @structexts=('',''); my @stdseqtexts=('',''); my @seqtexts=('',''); my @CMtexts=('','');
  
  
  for (my $i=0; $i<2; $i++) {
     while ($oristdseqs[$i]=~/\*\[\s*(\d+)\]\*/g) {
        my ($part_1_len, $part_2_len, $part_3_len, $spacing_len) = (length($`), $1, length($'), length($&));
        
        $stdseqtexts[$i].=$`.(' 'x$part_2_len);
        #reset orisedseqs here and after the loop, add back the final part of $oristdseqs...
        $oristdseqs[$i]=$';

        $oriseqs[$i]=~/(\S*?)\*\[\s*(\d+)\]\*/;
        $seqtexts[$i].=$1.(' 'x$part_2_len);
        $oriseqs[$i]=$';

        $oristrucs[$i]=~/(\S*?)\~+/;
        if ($part_1_len != length($1)) {die "error in .....\n";}
        $structexts[$i].=$1.(' 'x$part_2_len);
        $oristrucs[$i]=$';
        
        $oriCMs[$i]=~/(.{$part_1_len}).{$spacing_len}/; 
        $CMtexts[$i].=$1.(' 'x$part_2_len);
        $oriCMs[$i]=$';
     } 
        
     {
        $stdseqtexts[$i].=$oristdseqs[$i];
        $seqtexts[$i]   .=$oriseqs[$i];
        $structexts[$i] .=$oristrucs[$i];
        $CMtexts[$i]    .=$oriCMs[$i];
     }
     
     #add prefix spaces here
     for (my $j=$supposed_seq_start; $j<$best_hits[$i]->{'std_start'};$j++) {
        $stdseqtexts[$i]=' '.$stdseqtexts[$i];
        $structexts[$i] =' '.$structexts[$i];
        $CMtexts[$i]    =' '.$CMtexts[$i];
        $seqtexts[$i]   =' '.$seqtexts[$i];
     } 
        
     #add suffix spaces here
     for (my $k=$best_hits[$i]->{'std_end'}; $k<$supposed_seq_end; $k++) {
        $stdseqtexts[$i].=' ';
        $structexts[$i] .=' ';
        $CMtexts[$i]    .=' ';
        $seqtexts[$i]   .=' ';
     }

     $res->[$i]->{'stdseqtexts'} = $stdseqtexts[$i];
     $res->[$i]->{'structexts'}  = $structexts[$i];
     $res->[$i]->{'CMtexts'}     = $CMtexts[$i];
     $res->[$i]->{'seqtexts'}    = $seqtexts[$i];
  }

  return $res;
}

sub toString {
  my $self = shift;
  return join("\t", 'cm:', 
                    $self->cm1Bits,
                    $self->cm2Bits,
                    $self->alnLen,
                    $self->{'validStrucAln'},
                    "stem:",
                    $self->stLen,
                    $self->cmCount,
                    $self->pcmCount,
                    $self->stNoChange,
                    $self->stMismatch,
                    $self->stInfoCount,
                    "ss:",
                    $self->ssLen,
                    $self->ssNoChange,
                    $self->ssMismatch,
                    $self->ssInfoCount
             );
}

1;
