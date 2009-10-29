package EnsSliceTool;
use strict;
use Bio::Root::Root;
use Bio::PrimarySeq;
use vars qw($ID $VERSION @ISA @EXPORT);

@ISA=qw(Exporter);
@EXPORT=qw(&getSliceByGenomeDBRange
           &getStateByGadpGdbr 
           &getGenomeStateBySlice 
           &sliceToSeqObj 
           &getRepeatByGenomeDBRange
           &getWGAStateBySlice
           &getSliceAdpVersion
           &getRandomSlicesBySliceAdp
           &getHitContigStateBySlice
          );

$ID = 'EnsSliceTool';
$VERSION  = 0.1;

sub getRepeatByGenomeDBRange {
  my ($sliceadp, $dbrange) = (shift, shift);
  my $slice = getSliceByGenomeDBRange($sliceadp, $dbrange);
  if (!defined($slice)) {
     return undef;
  }
  
  return $slice->get_all_RepeatFeatures;
}

sub getSliceByGenomeDBRange {
  my ($sliceadp, $dbrange) = (shift, shift);
  if (!$dbrange->isa("genomedbrange")) {
     warn "$dbrange is not a genomedbrange\n";
     return undef;
  }

  my $newSlice = $sliceadp->fetch_by_region(
          $dbrange->getSeqLevel,
          $dbrange->getSeqRegion,
          $dbrange->getPos1,
          $dbrange->getPos2,
          $dbrange->getStrand,
         );

  if (!defined($newSlice)) {
     warn "Something wrong with ".$dbrange->getSequenceName."\n";
  }

  if ($newSlice->coord_system->version ne $dbrange->getDBVersion) {
     warn $dbrange->getDBVersion." is not compatible with the version in your sliceadaptor\n";
     return undef;
  } else {
     return $newSlice;
  }
}

sub getRandomSeqBySliceAdp {
  my $sliceadp = shift;
  my $length   = shift;
  
  
  return getRandomSeqsBySliceAdp($sliceadp, $length, 1);
}

sub getSliceVersion {
  my $slice = shift;
  if ($slice->isa("Bio::EnsEMBL::Slice")) {
     return $slice->coord_system->version;
  } else {
     warn "$slice is not a Bio::EnsEMBL::Slice object\n";
     return undef;
  }
}

sub getSliceAdpVersion {
  my $sliceadp = shift;
  if ($sliceadp->isa("Bio::EnsEMBL::DBSQL::SliceAdaptor")) {
     return $sliceadp->db->get_CoordSystemAdaptor()->fetch_by_name("chromosome")->version;
  } else {
     warn "$sliceadp is not a Bio::EnsEMBL::DBSQL::SliceAdaptor\n";
     return undef;
  }
}

sub getRandomSlicesBySliceAdp {
  my $sliceadp = shift;
  my $length   = shift;
  my $num_seqs = shift;

  my @seq_regs     = @{$sliceadp->fetch_all('toplevel')};
  my @new_seq_regs = ();
  foreach my $single_reg (@seq_regs) {
     if ($single_reg->seq_region_name eq "X" || $single_reg->seq_region_name eq "Y" || $single_reg->seq_region_name=~/^\d+$/) {
        push(@new_seq_regs, $single_reg);
     }
  }

  @seq_regs = @new_seq_regs;

  my %chr=();
  my %chr_accumulated_length = ();
  my $total_genome_length = 0;

  my @res = ();

  foreach my $single_region (@seq_regs) {
     #if ($single_region->seq_region_name eq "MT" || $single_region->seq_region_name eq "supercontig") {next;}
     $chr_accumulated_length{$single_region->seq_region_name} = $total_genome_length;
     $chr{$single_region->seq_region_name} = $single_region;
     $total_genome_length+=$single_region->length;
  }

  my $seq_count = 0;
  while ($seq_count < $num_seqs) {
     my $raw_pos = int(rand($total_genome_length));

     my ($chr, $pos) = _locate_chromosome_coord($raw_pos, 
                                                $length, 
                                                \%chr_accumulated_length, 
                                                \@seq_regs,
                                                \%chr);
     #print join("\t", $chr, $pos), "\n";
     if (defined($chr) && defined($pos)) {
         my $strand = (rand >= 0.5) ? 1 : -1;
         my $aSlice = $sliceadp->fetch_by_region(
                           "chromosome", 
                           $chr, 
                           $pos, 
                           $pos + $length - 1,
                           $strand);

         if ($aSlice->seq=~/N{4,}/) {next;}
         push(@res, $aSlice);
         $seq_count++;
     }
  }
  return \@res;
}

sub _locate_chromosome_coord {
  my $pos                        = shift;
  my $length                     = shift;
  my $chr_accumulated_length_ref = shift;
  my $seq_regs_ref               = shift;
  my $chr_ref                    = shift;

  my %chr_accumulated_length = %$chr_accumulated_length_ref;
  my @seq_regs               = @$seq_regs_ref;
  my %chr                    = %$chr_ref;

  my $in_chr = undef;
  my $offset = undef;

  foreach my $single_region (@seq_regs) {
     if ($pos > $chr_accumulated_length{$single_region->seq_region_name}) {
        $in_chr = $single_region->seq_region_name;
     } else {
        last;
     }
  }
  if (!defined($in_chr)) {return undef, undef;}
  $offset = $pos - $chr_accumulated_length{$in_chr};
  if (($offset + $length -1) > $chr{$in_chr}->length) {
     return undef, undef;
  }

  return $in_chr, $offset;
}

sub sliceToSeqObj {
  my $slice = shift;

  if (!$slice->isa("Bio::EnsEMBL::Slice")) {
     die "Your slice is not an object of Bio::EnsEMBL::Slice\n";
  }
  
  my $seqobj = new Bio::PrimarySeq(-display_id=>$slice->name, 
                                   -seq       =>$slice->seq);
  
  return $seqobj;
}

sub getGenomeStateByGadpGdbr {
  my ($sliceadp, $dbrange) = (shift, shift);
  my $newSlice = getSliceByGenomeDBRange($sliceadp, $dbrange);
  
  if (defined($newSlice)) {
     return getGenomeStateBySlice($newSlice);
  } else {
     warn "Something wrong with the genomedbrange you've provided\n";
     return undef;
  }
}

sub getWGAStateBySlice {
  my $slice = shift;
  my $res = "WGA";
  
  if (!$slice->isa("Bio::EnsEMBL::Slice")) {
     warn "Your slice is not an object of Bio::EnsEMBL::Slice";
     return "UN";
  }

  my $seqtext = $slice->seq();
  my $len = 0;
  if ($seqtext=~/(N{5,})/i) {
     $len = length($1);
  }

  if ($len > 0) {
     $res.="N$len";
  }

  return $res;
}

sub getHitContigStateBySlice {
  my $slice = shift;

  my @res;
  if (!$slice->isa("Bio::EnsEMBL::Slice")) {
     warn "Your slice is not an object of Bio::EnsEMBL::Slice";
     push(@res, "UN");
     return \@res;
  }
  
  my $project = $slice->project("contig");
  if (!defined($project)) {
     warn "Your slice cannot be projected to contigs\n";
     push(@res, "UN");
     return \@res;
  }

  foreach my $seg (@$project) {
     my $newSlice = $seg->to_Slice();
     my $seq_reg  = $newSlice->seq_region_name();
     my $name     = $seq_reg;

     my $hitSlice = $newSlice->adaptor->fetch_by_region("contig", $seq_reg);

     if ($seq_reg=~/(.*?)\./) {
        $name=$1;
     }

     my $res = IsWGA($name);

     {
        my $seqtext = $hitSlice->seq();
        my $len = 0;
        if ($seqtext=~/(N{5,})/i) {
           $len = length($1);
        }

        if ($len > 0) {
           $res.="N$len";
        } elsif ($res eq "CS") {
           $res = "FCS";
        }
     }
     push(@res, $res);
  }

  return \@res;
}

sub getGenomeStateBySlice {
  my $slice = shift;
  my $res;

  if (!$slice->isa("Bio::EnsEMBL::Slice")) {
     warn "Your slice is not an object of Bio::EnsEMBL::Slice";
     return "UN";
  }
  
  my $project = $slice->project("contig");
  if (!defined($project)) {
     warn "Your slice cannot be projected to contigs\n";
     return "UN";
  }
  
  foreach my $seg (@$project) {
     my $newSlice = $seg->to_Slice();
     my $seq_reg  = $newSlice->seq_region_name();
     my $name     = $seq_reg;

     if ($seq_reg=~/(.*?)\./) {
        $name=$1;
     }
     #print $slice->name, "\t", $name, "\n";

     $res = IsWGA($name);
     if ($res eq "WGA") {
        last;
     }
  }

  my $seqtext = $slice->seq();
  my $len = 0;
  if ($seqtext=~/(N{5,})/i) {
     $len = length($1);
  }

  if ($len > 0) {
     $res.="N$len";
  } elsif ($res eq "CS") {
     $res = "FCS";
  }
  
  return $res;
}

sub IsWGA_deprecated {
  my $name = shift;
  my $res  = `pfetch -D $name`;

  if ($res=~/no match/) {
     warn "Cannot find matches for $name\n";
     return "UN";
  }

  if ($res=~/shotgun/) {
     return "WGA";
  } else {
     return "CS";
  }
}

sub IsWGA {
  my $name = shift;
  if ($name=~/^CAA/i || $name=~/^AAD/i) {
     return "WGA";
  } else {
     return "CS";
  }
}

1;
