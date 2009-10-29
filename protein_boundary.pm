package protein_boundary;
use strict;
use Bio::Root::Root;
use Bio::PrimarySeq;
use gene_rec;
use simplerange;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'protein_boundary';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  #   hi gene   lo gene   hi gene   lo gene
  #   db1_gene1 db1_gene2 db2_gene1 db2_gene2 
  my ($ori_slice, $gene1_1, $gene1_2, $gene2_1, $gene2_2, $strand, $ensDB1, $ensDB2) = 
     $self->_rearrange([qw(ORI_SLICE GENE1_1 GENE1_2 GENE2_1 GENE2_2 STRAND DB1 DB2)], @args);

  if ($gene1_1->isa('gene_rec') && $gene1_2->isa('gene_rec') && 
      $gene2_1->isa('gene_rec') && $gene2_2->isa('gene_rec')) {

      $self->{'gene1_1'} = $gene1_1;
      $self->{'gene1_2'} = $gene1_2;
      $self->{'gene2_1'} = $gene2_1;
      $self->{'gene2_2'} = $gene2_2;

      $self->{'ori_slice'} = $ori_slice;

      my $DB1region = $self->DB1BoundedRegion();
      my $DB2region = $self->DB2BoundedRegion();

      if ((!defined($DB1region)) || (!defined($DB2region))) {
         warn "not on the same chromosome\n";
         return undef;
      }

      $self->{'slice1'}  = 
         $ensDB1->get_SliceAdaptor()->fetch_by_region('chromosome', 
                                                      $DB1region->getSeqName,
                                                      $DB1region->getPos1,
                                                      $DB1region->getPos2
                                                     );
      $self->{'slice2'}  = 
         $ensDB2->get_SliceAdaptor()->fetch_by_region('chromosome',
                                                      $DB2region->getSeqName, 
                                                      $DB2region->getPos1,
                                                      $DB2region->getPos2
                                                     );

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
  
  my $intergenic_chr   = $hi_gene->chr_name;
  if ($intergenic_chr ne $lo_gene->chr_name) {
     return undef;
  }

  my $intergenic_end   = $hi_gene->end;
  my $intergenic_start = $lo_gene->start;

  return new simplerange(-seqname=>$intergenic_chr, -pos1=>$intergenic_start, -pos2=>$intergenic_end);
}

sub DB2BoundedRegion {
  my $self = shift;
  my $hi_gene = $self->DB2HiGene;
  my $lo_gene = $self->DB2LoGene;
  
  my $intergenic_chr   = $hi_gene->chr_name;
  if ($intergenic_chr ne $lo_gene->chr_name) {
     warn "On different chromosome\n";
     return undef;
  }

  my $intergenic_end   = $hi_gene->end;
  my $intergenic_start = $lo_gene->start;

  return new simplerange(-seqname=>$intergenic_chr, -pos1=>$intergenic_start, -pos2=>$intergenic_end);
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
