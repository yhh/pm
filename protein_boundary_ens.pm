package protein_boundary_ens;
use strict;
use Bio::Root::Root;
use Bio::PrimarySeq;
use simplerange;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'protein_boundary_ens';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  #   hi gene   lo gene   hi gene   lo gene
  #   5'	3'	  5'	    3'
  #   db1_gene1 db1_gene2 db2_gene1 db2_gene2 
  my ($ori_slice, $gene1_1, $gene1_2, $gene2_1, $gene2_2, $strand, $isComplicated, $db1desc, $db2desc) = 
  $self->_rearrange([qw(ORI_SLICE GENE1_1 GENE1_2 GENE2_1 GENE2_2 STRAND COMPLICATED DB1DESC DB2DESC)], @args);

  if ($gene1_1->isa('Bio::EnsEMBL::Gene') && $gene1_2->isa('Bio::EnsEMBL::Gene') && 
      $gene2_1->isa('Bio::EnsEMBL::Gene') && $gene2_2->isa('Bio::EnsEMBL::Gene')) {

      $self->{'gene1_1'} = $gene1_1;
      $self->{'gene1_2'} = $gene1_2;
      $self->{'gene2_1'} = $gene2_1;
      $self->{'gene2_2'} = $gene2_2;

      $self->{'db1desc'} = $db1desc;
      $self->{'db2desc'} = $db2desc; 

      $self->{'ori_slice'} = $ori_slice;

      $self->{'complicated'} = $isComplicated;

      my $DB1region = $self->DB1BoundedRegion();
      my $DB2region = $self->DB2BoundedRegion();

      my $ensDB1 = $gene1_1->slice->adaptor->db;
      my $ensDB2 = $gene2_1->slice->adaptor->db;

      #if ((!defined($DB1region)) || (!defined($DB2region))) {
      #   warn "not on the same chromosome\n";
      #   return undef;
      #}

      if (defined($DB1region)) {
         $self->{'slice1'}  = 
         $ensDB1->get_SliceAdaptor()->fetch_by_region('chromosome', 
                                                      $DB1region->getSeqName,
                                                      $DB1region->getPos1,
                                                      $DB1region->getPos2
                                                     );
         $self->{'isDB1SameChr'} = 1;
      } else {
         $self->{'isDB1SameChr'} = 0;
      }

      if (defined($DB2region)) {
         $self->{'slice2'}  = 
         $ensDB2->get_SliceAdaptor()->fetch_by_region('chromosome',
                                                      $DB2region->getSeqName, 
                                                      $DB2region->getPos1,
                                                      $DB2region->getPos2
                                                     );
         $self->{'isDB2SameChr'} = 1;
      } else {
         $self->{'isDB2SameChr'} = 0;
      }

      if (!defined($strand)) {
         my $strand_check_1 = $gene1_1->strand * $gene2_1->strand;
         my $strand_check_2 = $gene1_2->strand * $gene2_2->strand;
         if ($strand_check_1 == 0) {
            $strand = 0;
         } elsif ($strand_check_1 == $strand_check_2) {
            if ($gene1_1->strand != $gene2_1->strand) {
               $strand = -1;
            } else {
               $strand = 1;
            }
         } else {
            $strand = -5;
         }
      }
      $self->{'strand'}  = $strand;
  } else {
     die 'gene1_1 or gene1_2 or gene2_1 or gene2_2 is not a gene_rec'."\n";
     return undef;
  }

  return $self;
}

sub toString {
  my $self = shift;

  return     join("\t", $self->DB1Slice->name, 
                        $self->getDB1Desc,
			$self->getStringlizedDB1BoundedRegionOverlap,
			$self->stringlizedIsComplicated,
                        $self->DB2Slice->name,
                        $self->getDB2Desc,
			$self->getStrand,
			$self->getStringlizedDB2BoundedRegionOverlap
                 );
}

sub getDB1Desc {
  my $self = shift;
  my $currentDesc = $self->{'db1desc'};
  if (!defined($currentDesc) || $currentDesc eq "") {
     return "normal";
  } else {
     return $currentDesc;
  }
}

sub getDB2Desc {
  my $self = shift;
  my $currentDesc = $self->{'db2desc'}; 
  if (!defined($currentDesc) || $currentDesc eq "") {
     return "normal";
  } else {
     return $currentDesc;
  }
}

sub setDB1Desc {
  my $self = shift;
  my $new_desc = shift;
  $self->{'db1desc'} = $new_desc;
}

sub setDB2Desc {
  my $self = shift;
  my $new_desc = shift;
  $self->{'db2desc'} = $new_desc;
}

sub isDB1SameChr {
  my $self = shift;
  return $self->{'isDB1SameChr'};
}

sub isDB2SameChr {
  my $self = shift;
  return $self->{'isDB2SameChr'};
}

sub isComplicated {
  my $self = shift;
  return $self->{'complicated'};
}

sub stringlizedIsComplicated {
  my $self = shift;
  return ($self->isComplicated ? "complicated" : "simple");
}

sub getDB1BoundedRegionOverlap {
  my $self = shift;
  return $self->{'DB1BoundedRegionOverlap'};
}

sub getStringlizedDB1BoundedRegionOverlap {
  my $self = shift;
  return (defined($self->getDB1BoundedRegionOverlap) ? "DB1BoundaryOverlap" : "DB1BoundaryOK");
}

sub getDB2BoundedRegionOverlap {
  my $self = shift;
  return $self->{'DB2BoundedRegionOverlap'};
}

sub getStringlizedDB2BoundedRegionOverlap {
  my $self = shift;
  return (defined($self->getDB2BoundedRegionOverlap) ? "DB2BoundaryOverlap" : "DB2BoundaryOK");
}

sub ori_slice {
  my $self = shift;
  return $self->{'ori_slice'};
}

sub DB1Slice {
  my $self = shift;
  return $self->{'slice1'};
}

sub DB2Slice {
  my $self = shift;
  return $self->{'slice2'};
}

sub ori_slice_seq {
  my $self = shift;
  my $ori_slice = $self->ori_slice;
  return Bio::PrimarySeq->new(-seq=>$ori_slice->seq(), -display_id=>$ori_slice->name());
}

sub DB1Seq {
  my $self = shift;
  my $db1slice = $self->DB1Slice;
  return Bio::PrimarySeq->new(-seq=>$db1slice->seq(), -display_id=>$db1slice->name());
}

sub DB2Seq {
  my $self = shift;
  my $db2slice = $self->DB2Slice;
  return Bio::PrimarySeq->new(-seq=>$db2slice->seq(), -display_id=>$db2slice->name());
}

sub DB1GenePair {
  my $self = shift;
  my $pair = {-gene1=>$self->{'gene1_1'}, -gene2=>$self->{'gene1_2'}};
  return $pair;
}

sub DB2GenePair {
  my $self = shift;
  my $pair = {-gene1=>$self->{'gene2_1'}, -gene2=>$self->{'gene2_2'}};
  return $pair;
}

sub DB1HiGene {
  my $self = shift;
  return $self->{'gene1_1'};
}

sub DB1LoGene {
  my $self = shift;
  return $self->{'gene1_2'};
}

sub DB1BoundedRegion {
  my $self = shift;
  my $hi_gene = $self->DB1HiGene;
  my $lo_gene = $self->DB1LoGene;
  
  my $intergenic_chr   = $hi_gene->slice->seq_region_name;
  if ($intergenic_chr ne $lo_gene->slice->seq_region_name) {
     $self->{'DB1BoundedRegionOverlap'} = undef;
     return undef;
  }

  my $intergenic_start = $hi_gene->feature_Slice->end;
  my $intergenic_end   = $lo_gene->feature_Slice->start;

  my $simplerange1 = new simplerange(-seqname=>$intergenic_chr, -pos1=>$hi_gene->feature_Slice->start, 
                                                                -pos2=>$hi_gene->feature_Slice->end);
  my $simplerange2 = new simplerange(-seqname=>$intergenic_chr, -pos1=>$lo_gene->feature_Slice->start,
                                                                -pos2=>$lo_gene->feature_Slice->start);

  $self->{'DB1BoundedRegionOverlap'} = $simplerange1->checkOverlap($simplerange2);

  #if (!defined($self->{'DB1BoundedRegionOverlap'})) {
     #return new simplerange(-seqname=>$intergenic_chr, -pos1=>$intergenic_start, -pos2=>$intergenic_end);

  #for coping genes like snoRNAs and miRNAs which may be in introns....

     return new simplerange(-seqname=>$intergenic_chr, -pos1=>$hi_gene->feature_Slice->start, 
                                                       -pos2=>$lo_gene->feature_Slice->end);
  #} else {
  #   return new simplerange(-seqname=>$intergenic_chr, -pos1=>$hi_gene->feature_Slice->start, 
  #                                                     -pos2=>$lo_gene->feature_Slice->end);
  #}
}

sub DB2BoundedRegion {
  my $self = shift;
  my $hi_gene = $self->DB2HiGene;
  my $lo_gene = $self->DB2LoGene;
  
  my $intergenic_chr   = $hi_gene->slice->seq_region_name;
  if ($intergenic_chr ne $lo_gene->slice->seq_region_name) {
     $self->{'DB2BoundedRegionOverlap'} = undef;
     return undef;
  }

  my $intergenic_start = $hi_gene->feature_Slice->end;
  my $intergenic_end   = $lo_gene->feature_Slice->start;


  my $simplerange1 = new simplerange(-seqname=>$intergenic_chr, -pos1=>$hi_gene->feature_Slice->start, 
                                                                -pos2=>$hi_gene->feature_Slice->end);
  my $simplerange2 = new simplerange(-seqname=>$intergenic_chr, -pos1=>$lo_gene->feature_Slice->start,
                                                                -pos2=>$lo_gene->feature_Slice->start);

  $self->{'DB2BoundedRegionOverlap'} = $simplerange1->checkOverlap($simplerange2);

  #if (!defined($self->{'DB2BoundedRegionOverlap'})) {
  #   return new simplerange(-seqname=>$intergenic_chr, -pos1=>$intergenic_start, -pos2=>$intergenic_end);
  #} else {
  
  #for coping genes like snoRNAs and miRNAs which may be in introns....
  
     return new simplerange(-seqname=>$intergenic_chr, -pos1=>$hi_gene->feature_Slice->start, 
                                                       -pos2=>$lo_gene->feature_Slice->end);
  #}
}

sub DB2HiGene {
  my $self = shift;
  return $self->{'gene2_1'};
}

sub DB2LoGene {
  my $self = shift;
  return $self->{'gene2_2'};
}

sub getStrand {
  my $self = shift;
  return $self->{'strand'};
}

1;
