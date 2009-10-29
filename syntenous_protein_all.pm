package syntenous_protein_all;
use strict;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::Root::Root;
use sortfuncs;
use protein_boundary_ens;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'syntenous_protein_all';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($conf_file, $host, $dbuser, $port, $dbname, $assembly1, $assembly2, $species1, $species2) = 
  $self->_rearrange([qw(CONF_FILE HOST DBUSER PORT DBNAME ASSEMBLY1 ASSEMBLY2 SPECIES1 SPECIES2)], @args);
  
  if (!defined($port)) {$port = 3306;}
  my $db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor (-host => $host,
						        -user => $dbuser,
						        -port => $port,
						        -dbname => $dbname,
						        -conf_file => $conf_file);
  
  my @dbconf = (undef, undef);
  my $server_conf = do $conf_file || die "Cannot find $conf_file or config file err!\n";

  foreach my $single_species (@$server_conf) {
     my ($species_name, $species_version, $species_conf) = @$single_species;
     #print %$species_conf, "\n";
     if ($species_name eq $species1 && $species_version eq $assembly1) {
        $dbconf[0] = $species_conf;
     } elsif ($species_name eq $species2 && $species_version eq $assembly2) {
        $dbconf[1] = $species_conf;
     }
  }

  my @ensdb = (undef, undef);
  for (my $i=0; $i<2; $i++) {
     if (!defined($dbconf[$i]->{'port'})) {$dbconf[$i]->{'port'} = 3306;}
     $ensdb[$i] = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbconf[$i]->{'host'},
                                                     -user   => $dbconf[$i]->{'user'},
                                                     -dbname => $dbconf[$i]->{'dbname'},
                                                     -port   => $dbconf[$i]->{'port'}
                                                    );
  }

  $self->{'assembly1'} = $assembly1;
  $self->{'assembly2'} = $assembly2;

  $self->{'species1'} = $species1;
  $self->{'species2'} = $species2;

  $self->{'db'}     = $db;
  $self->{'ensdb1'} = $ensdb[0];
  $self->{'ensdb2'} = $ensdb[1];

  $self->{'genead1'}= $ensdb[0]->get_GeneAdaptor();
  $self->{'genead2'}= $ensdb[1]->get_GeneAdaptor();

  $self->{'sliceadaptor1'} = $db->get_db_adaptor($species1,$assembly1)->get_SliceAdaptor;
  $self->{'sliceadaptor2'} = $db->get_db_adaptor($species2,$assembly2)->get_SliceAdaptor;
 
  warn "pre getting genomedb\n";
  my $gdadp        = $db->get_GenomeDBAdaptor;
  my $gdb1         = $gdadp->fetch_by_name_assembly($species1, $assembly1);
  my $gdb2         = $gdadp->fetch_by_name_assembly($species2, $assembly2);

  $self->{'genome_db_id1'} = $gdb1->dbID;
  $self->{'genome_db_id2'} = $gdb2->dbID;
  $self->{'homoloadp'}     = $db->get_HomologyAdaptor;

  warn "pre fetching top level\n";
  my @seq_regs= @{$self->{'sliceadaptor1'}->fetch_all('toplevel')};

  my %chr=();
  foreach my $seq_r(@seq_regs) {
    $chr{$seq_r->seq_region_name} = $seq_r;
  }

  $self->{'chr_regs'} = \%chr;

  warn "pre fetching all by genome pairs\n";
  my $homologyListRef = $self->{'homoloadp'}->fetch_all_by_genome_pair($self->{'genome_db_id1'}, 
                                                                       $self->{'genome_db_id2'});
  
  warn "post fetching all by genome pairs\n";
  warn "pre fetching UBRHs\n";
  my ($homolog_pairs_ref1, $homolog_pairs_ref2) = $self->getUBRHs($homologyListRef);

  warn "post fetching UBRHs\n";
  my @ordered_raw_genome_db_ranges_by_db1 = sort orderbygenomedbrange_raw keys %$homolog_pairs_ref1;
  my @ordered_raw_genome_db_ranges_by_db2 = sort orderbygenomedbrange_raw keys %$homolog_pairs_ref2;

  warn "post ordering dbrange1 and dbrange2\n";

  $self->{'homolog_pairs_ref1'} = $homolog_pairs_ref1;
  $self->{'homolog_pairs_ref2'} = $homolog_pairs_ref2;
  
  $self->{'ordered_raw_genome_db_ranges_by_db1'} = \@ordered_raw_genome_db_ranges_by_db1;
  $self->{'ordered_raw_genome_db_ranges_by_db2'} = \@ordered_raw_genome_db_ranges_by_db2;

  return $self; # success - we hope!
}

sub getBoundaryByDB1 {
  my $self= shift;
  my $genomedbrange = shift;

  # hack for skipping the intervening proteins from chromosomes 
  # different to 5' and 3' flanking syntenous proteins
  
  my $ignoreIntervenSingleton = shift;

  if (!defined($ignoreIntervenSingleton)) {
     $ignoreIntervenSingleton = 1;
  }

  if (!$genomedbrange->isa("genomedbrange")) {return undef, undef;}
  my $ordered_raw_genome_db_ranges_by_db1_ref = $self->{'ordered_raw_genome_db_ranges_by_db1'};

  my $loweridx  = dividefinder(0,              scalar(@$ordered_raw_genome_db_ranges_by_db1_ref)-1, 
                              $genomedbrange, $ordered_raw_genome_db_ranges_by_db1_ref);
  
  my ($pLowerIdx, $pUpperIdx) = $self->getDB1BoundariesIdx($genomedbrange, $loweridx, $ignoreIntervenSingleton);
  
  my $ori_slice = $self->{'sliceadaptor1'}->fetch_by_region($genomedbrange->getSeqLevel, 
                                                            $genomedbrange->getSeqRegion, 
                                                            $genomedbrange->getPos1, 
                                                            $genomedbrange->getPos2, 
                                                            $genomedbrange->getStrand);

  my @hiGenePair = @{$self->{'homolog_pairs_ref1'}->{$ordered_raw_genome_db_ranges_by_db1_ref->[$pLowerIdx]}};
  my @loGenePair = @{$self->{'homolog_pairs_ref1'}->{$ordered_raw_genome_db_ranges_by_db1_ref->[$pUpperIdx]}};

  my $newBoundary = new protein_boundary_ens(-ORI_SLICE   => $ori_slice,
                                             -GENE1_1     => $hiGenePair[0], 
                                             -GENE2_1     => $hiGenePair[1],
                                             -GENE1_2     => $loGenePair[0],
                                             -GENE2_2     => $loGenePair[1],
                                             -COMPLICATED => (abs($pUpperIdx - $pLowerIdx) > 1)
                                            );
  return $newBoundary;
}

sub getBoundaryByDB1_ext {
  my $self= shift;
  my $genomedbrange = shift;

  # hack for skipping the intervening proteins from chromosomes 
  # different to 5' and 3' flanking syntenous proteins
  
  my $ignoreIntervenSingleton = shift;

  if (!defined($ignoreIntervenSingleton)) {
     $ignoreIntervenSingleton = 1;
  }

  if (!$genomedbrange->isa("genomedbrange")) {return undef, undef;}
  my $ordered_raw_genome_db_ranges_by_db1_ref = $self->{'ordered_raw_genome_db_ranges_by_db1'};

  my $loweridx  = dividefinder(0,              scalar(@$ordered_raw_genome_db_ranges_by_db1_ref)-1, 
                              $genomedbrange, $ordered_raw_genome_db_ranges_by_db1_ref);
  
  my ($pLowerIdx, $pUpperIdx) = $self->getDB1BoundariesIdx($genomedbrange, $loweridx, $ignoreIntervenSingleton);
  
  #my $ori_slice = $self->{'sliceadaptor1'}->fetch_by_region($genomedbrange->getSeqLevel, 
  #                                                          $genomedbrange->getSeqRegion, 
  #                                                          $genomedbrange->getPos1, 
  #                                                          $genomedbrange->getPos2, 
  #                                                          $genomedbrange->getStrand);

  my @hiGenePair;
  if ($pLowerIdx == -5) {
     @hiGenePair = (undef, undef);
  } else {
     @hiGenePair = @{$self->{'homolog_pairs_ref1'}->{$ordered_raw_genome_db_ranges_by_db1_ref->[$pLowerIdx]}};
  }

  my @loGenePair;
  if ($pUpperIdx == -5) {
     @loGenePair = (undef, undef);
  } else {
     @loGenePair = @{$self->{'homolog_pairs_ref1'}->{$ordered_raw_genome_db_ranges_by_db1_ref->[$pUpperIdx]}};
  }

  return @hiGenePair, @loGenePair, (abs($pUpperIdx - $pLowerIdx) > 1);
}

sub getBoundaryByDB2 {
  my $self= shift;
  my $genomedbrange = shift;
  my $ignoreIntervenSingleton = shift;

  if (!defined($ignoreIntervenSingleton)) {
     $ignoreIntervenSingleton = 1;
  }

  if (!$genomedbrange->isa("genomedbrange")) {return undef, undef;}
  my $ordered_raw_genome_db_ranges_by_db2_ref = $self->{'ordered_raw_genome_db_ranges_by_db2'};

  my $loweridx  = dividefinder(0,              scalar(@$ordered_raw_genome_db_ranges_by_db2_ref)-1, 
                              $genomedbrange, $ordered_raw_genome_db_ranges_by_db2_ref);

  my ($pLowerIdx, $pUpperIdx) = $self->getDB2BoundariesIdx($loweridx, $ignoreIntervenSingleton);

  my $ori_slice = $self->{'sliceadaptor1'}->fetch_by_region($genomedbrange->getSeqLevel, 
                                                            $genomedbrange->getSeqRegion, 
                                                            $genomedbrange->getPos1, 
                                                            $genomedbrange->getPos2, 
                                                            $genomedbrange->getStrand);

  my @hiGenePair = @{$self->{'homolog_pairs_ref2'}->{$ordered_raw_genome_db_ranges_by_db2_ref->[$pLowerIdx]}};
  my @loGenePair = @{$self->{'homolog_pairs_ref2'}->{$ordered_raw_genome_db_ranges_by_db2_ref->[$pUpperIdx]}};

  my $newBoundary = new protein_boundary_ens(-ORI_SLICE   => $ori_slice,
                                             -GENE1_1     => $hiGenePair[0], 
                                             -GENE2_1     => $hiGenePair[1],
                                             -GENE1_2     => $loGenePair[0],
                                             -GENE2_2     => $loGenePair[1],
                                             -COMPLICATED => (abs($pUpperIdx - $pLowerIdx) > 1)
                                            );
  return $newBoundary;
  
}

sub dividefinder {
  #this is a private function which is not called outside this package
  my ($lower, $upper, $genomedbrange, $genomedbrangesRef) = (shift, shift, shift, shift);

  if ($upper - $lower == 1) {
     if ($genomedbrange->compareRaw($genomedbrangesRef->[$upper]) > 0) {
        return $upper+1;
     } elsif ($genomedbrange->compareRaw($genomedbrangesRef->[$lower]) < 0) {
        return $lower-1;
     } else {
        return $lower;
     }
  }

  my $middle=int(($upper+$lower)/2);
  #print join("\t", $upper, $middle, $lower), " upper middle lower\n";
  #if (($genomedbrange->compareRaw($genomedbrangesRef->[$lower])  >=0) && 
  #    ($genomedbrange->compareRaw($genomedbrangesRef->[$middle]) < 0)) {
  if ($genomedbrange->compareRaw($genomedbrangesRef->[$middle]) < 0) {

     return dividefinder($lower, $middle, $genomedbrange, $genomedbrangesRef);
  } else {
     return dividefinder($middle, $upper, $genomedbrange, $genomedbrangesRef);
  }
}

sub getDB1BoundariesIdx {
  my $self = shift;
  my $ori_gdbr = shift;
  my $lowerIdx = shift;
  my $ignoreInterven = shift;
  
  if (!defined($ignoreInterven)) {$ignoreInterven = 0;}
  if (!$ignoreInterven) {return $lowerIdx, $lowerIdx+1;}

  my ($newLowerIdx, $newUpperIdx) = ($lowerIdx, $lowerIdx + 1);

  my $ordered_raw_genome_db_ranges_by_db1_ref = $self->{'ordered_raw_genome_db_ranges_by_db1'};

  # unmarked by Yen-Hua Huang, 05/10/06
  my $ori_db1_chr = $ori_gdbr->getSeqRegion();
  
  while ($newLowerIdx >= 0) {
     my $prev_raw_genomedb1range = $ordered_raw_genome_db_ranges_by_db1_ref->[$newLowerIdx - 1];
     my $this_raw_genomedb1range = $ordered_raw_genome_db_ranges_by_db1_ref->[$newLowerIdx];
     my $next_raw_genomedb1range = $ordered_raw_genome_db_ranges_by_db1_ref->[$newLowerIdx + 1];

     if ($ori_db1_chr ne $self->{'homolog_pairs_ref1'}->{$prev_raw_genomedb1range}->[0]->seq_region_name) {
        #return -5, $newlowerIdx + 1;
        $newLowerIdx = -5;
        last;
     }

     my ($prev_seqregion, $this_seqregion, $next_seqregion) = 
             ($self->{'homolog_pairs_ref1'}->{$prev_raw_genomedb1range}->[1]->seq_region_name, 
              $self->{'homolog_pairs_ref1'}->{$this_raw_genomedb1range}->[1]->seq_region_name,
              $self->{'homolog_pairs_ref1'}->{$next_raw_genomedb1range}->[1]->seq_region_name);
     
     my $lowerOK = 
        #($self->{'homolog_pairs_ref1'}->{$this_raw_genomedb1range}->[1]->coord_system_name eq "chromosome" && 
        ($this_seqregion ne "" && 
         ($this_seqregion eq $prev_seqregion || $this_seqregion eq $next_seqregion)
        ); 

     if (defined($ori_gdbr->checkOverlapRaw($this_raw_genomedb1range))) {
        $lowerOK = 0;
     }

     if ($lowerOK) {last;}
     $newLowerIdx--;
  }
 
  while ($newUpperIdx < scalar(@$ordered_raw_genome_db_ranges_by_db1_ref)) {
     my $prev_raw_genomedb1range = $ordered_raw_genome_db_ranges_by_db1_ref->[$newUpperIdx - 1];
     my $this_raw_genomedb1range = $ordered_raw_genome_db_ranges_by_db1_ref->[$newUpperIdx];
     my $next_raw_genomedb1range = $ordered_raw_genome_db_ranges_by_db1_ref->[$newUpperIdx + 1];

     # new modification here by Yen-Hua 05/10/06
     if ($ori_db1_chr ne $self->{'homolog_pairs_ref1'}->{$next_raw_genomedb1range}->[0]->seq_region_name) {
        #return $newLowerIdx, -5;
        $newUpperIdx = -5;
        last;
     }

     my ($prev_seqregion, $this_seqregion, $next_seqregion) = 
             ($self->{'homolog_pairs_ref1'}->{$prev_raw_genomedb1range}->[1]->seq_region_name, 
              $self->{'homolog_pairs_ref1'}->{$this_raw_genomedb1range}->[1]->seq_region_name,
              $self->{'homolog_pairs_ref1'}->{$next_raw_genomedb1range}->[1]->seq_region_name);
     
     my $UpperOK = 
        #($self->{'homolog_pairs_ref1'}->{$this_raw_genomedb1range}->[1]->coord_system_name eq "chromosome" &&
        ($this_seqregion ne "" && 
         ($this_seqregion eq $prev_seqregion || $this_seqregion eq $next_seqregion)
        ); 

     if (defined($ori_gdbr->checkOverlapRaw($this_raw_genomedb1range))) {
        $UpperOK = 0;
     }

     if ($UpperOK) {last;}
     $newUpperIdx++;
  }
  return $newLowerIdx, $newUpperIdx;
}

sub getDB2BoundariesIdx {
  my $self = shift;
  my $lowerIdx = shift;
  my $ignoreInterven = shift;
  
  if (!defined($ignoreInterven)) {$ignoreInterven = 0;}
  if (!$ignoreInterven) {return $lowerIdx, $lowerIdx+1;}

  my ($newLowerIdx, $newUpperIdx) = ($lowerIdx, $lowerIdx + 1);

  my $ordered_raw_genome_db_ranges_by_db2_ref = $self->{'ordered_raw_genome_db_ranges_by_db2'};
  
  while ($newLowerIdx >= 0) {
     my $prev_raw_genomedb2range = $ordered_raw_genome_db_ranges_by_db2_ref->[$newLowerIdx - 1];
     my $this_raw_genomedb2range = $ordered_raw_genome_db_ranges_by_db2_ref->[$newLowerIdx];
     my $next_raw_genomedb2range = $ordered_raw_genome_db_ranges_by_db2_ref->[$newLowerIdx + 1];

     my ($prev_seqregion, $this_seqregion, $next_seqregion) = 
             ($self->{'homolog_pairs_ref2'}->{$prev_raw_genomedb2range}->[0]->seq_region_name, 
              $self->{'homolog_pairs_ref2'}->{$this_raw_genomedb2range}->[0]->seq_region_name,
              $self->{'homolog_pairs_ref2'}->{$next_raw_genomedb2range}->[0]->seq_region_name);
     
     my $lowerOK = ($this_seqregion ne "" && ($this_seqregion eq $prev_seqregion || 
                                              $this_seqregion eq $next_seqregion)
                   ); 

     if ($lowerOK) {last;}
     $newLowerIdx--;
  }
 
  while ($newUpperIdx < scalar(@$ordered_raw_genome_db_ranges_by_db2_ref)) {
     my $prev_raw_genomedb2range = $ordered_raw_genome_db_ranges_by_db2_ref->[$newUpperIdx - 1];
     my $this_raw_genomedb2range = $ordered_raw_genome_db_ranges_by_db2_ref->[$newUpperIdx];
     my $next_raw_genomedb2range = $ordered_raw_genome_db_ranges_by_db2_ref->[$newUpperIdx + 1];

     my ($prev_seqregion, $this_seqregion, $next_seqregion) = 
             ($self->{'homolog_pairs_ref2'}->{$prev_raw_genomedb2range}->[0]->seq_region_name, 
              $self->{'homolog_pairs_ref2'}->{$this_raw_genomedb2range}->[0]->seq_region_name,
              $self->{'homolog_pairs_ref2'}->{$next_raw_genomedb2range}->[0]->seq_region_name);
     
     my $UpperOK = ($this_seqregion ne "" && ($this_seqregion eq $prev_seqregion || 
                                              $this_seqregion eq $next_seqregion)
                   ); 

     if ($UpperOK) {last;}
     $newUpperIdx++;
  }
  return $newLowerIdx, $newUpperIdx;
}

sub getAllHomologPairsOrderedByDB1Range {
  my $self = shift;
  my @res = ();
  foreach my $single_db_range (@{$self->{'ordered_raw_genome_db_ranges_by_db1'}}) {
     push(@res, $self->{'homolog_pairs_ref1'}->{$single_db_range});
  }
  return \@res;
}

sub getAllHomologPairsOrderedByDB2Range {
  my $self = shift;
  my @res = ();
  foreach my $single_db_range (@{$self->{'ordered_raw_genome_db_ranges_by_db2'}}) {
     push(@res, $self->{'homolog_pairs_ref2'}->{$single_db_range});
  }
  return \@res;
}

sub getDB1HomologGenes {
  my $self = shift;
  my @genes = ();
  foreach my $single_db_range (@{$self->{'ordered_raw_genome_db_ranges_by_db1'}}) {
     push(@genes, $self->{'homolog_pairs_ref1'}->{$single_db_range}->[0]);
  }
  return \@genes;
}

sub getDB2HomologGenes {
  my $self = shift;
  my @genes = ();
  foreach my $single_db_range (@{$self->{'ordered_raw_genome_db_ranges_by_db2'}}) {
     push(@genes, $self->{'homolog_pairs_ref1'}->{$single_db_range}->[1]);
  }
  return \@genes;
}

sub ensDB1 {
  my $self = shift;
  return $self->{'ensdb1'};
}

sub ensDB2 {
  my $self = shift;
  return $self->{'ensdb2'};
}

sub getUBRHs {
  my $self = shift;
  my $homology_list_ref = shift;
  my $returnHashRef1 = {};
  my $returnHashRef2 = {};
  foreach my $single_homolog (@$homology_list_ref) {
	if ($single_homolog->description ne "UBRH") {next;}
	my @homologmemattrs = @{$single_homolog->get_Member_Attribute_by_source("ENSEMBLGENE")};
	my $homologcount = scalar(@homologmemattrs);
	if ($homologcount ==2 ) {
		my @members          = ($homologmemattrs[0]->[0], $homologmemattrs[1]->[0]);
		my $members_hash_ref = $self->sortByDBOrder(\@members);
		my @syntenous_genes  = ($members_hash_ref->{$self->{'genome_db_id1'}}->get_Gene, 
                                        $members_hash_ref->{$self->{'genome_db_id2'}}->get_Gene);

		#hack for genes on supercontigs in DB2
		if ($syntenous_genes[1]->feature_Slice->coord_system_name ne "chromosome") {
                   next;
                }

		$returnHashRef1->{$syntenous_genes[0]->feature_Slice->name} = \@syntenous_genes;
		$returnHashRef2->{$syntenous_genes[1]->feature_Slice->name} = \@syntenous_genes;
	} 
	else {
		warn $homologcount, "-----------\n";
	}
  }
  return $returnHashRef1, $returnHashRef2;
}

sub sortByDBOrder {
  my $self = shift;
  my $members_ref  = shift;
  my %members = ();
  foreach my $single_member (@$members_ref) {
     if ($single_member->genome_db_id eq $self->{'genome_db_id1'}) {
     	$members{$self->{'genome_db_id1'}} = $single_member;
     } elsif ($single_member->genome_db_id eq $self->{'genome_db_id2'}) {
	$members{$self->{'genome_db_id2'}} = $single_member;
     }
  }
  return \%members;
}

#sub extendDB1ChrEndBoundary {
#  my $self = shift;
#  my $target_genomedbrange = shift;
#  
#}

sub extendDB1Boundary {
  my $self = shift;
  my $oldBoundary = shift;
  my $target_genomedbrange = shift;
  
  my $homolog_pairs_ref1 = $self->{'homolog_pairs_ref1'};

  my $lower_genomedb1range = new genomedbrange(-dbrange=>$oldBoundary->DB1HiGene->feature_Slice->name);

  my $lower_genomedb2range = new genomedbrange(
     -dbrange=>$homolog_pairs_ref1->{$lower_genomedb1range->getSequenceName}->[1]->feature_Slice->name);

  my $ordered_raw_genome_db_ranges_by_db1_ref = $self->{'ordered_raw_genome_db_ranges_by_db1'};

  my $loweridx  = dividefinder(0,                     scalar(@$ordered_raw_genome_db_ranges_by_db1_ref)-1, 
                               $lower_genomedb1range, $ordered_raw_genome_db_ranges_by_db1_ref);
  
  my $guessed_direction_for_extension = guessForwardDB2Direction(
                              $loweridx, 
                              $ordered_raw_genome_db_ranges_by_db1_ref,
                              $homolog_pairs_ref1
                             );

  my $guessed_distance = guessDistance($lower_genomedb1range, $target_genomedbrange);
  my $new_db2range      = $lower_genomedb2range->shiftDistance($guessed_distance, $guessed_direction_for_extension);
  my $new_db1range      = $lower_genomedb1range->shiftDistance($guessed_distance, 1);


  my $pseudo_sliceDB2 = $self->{'sliceadaptor2'}->fetch_by_region($new_db2range->getSeqLevel, 
                                                                  $new_db2range->getSeqRegion, 
                                                                  $new_db2range->getPos1, 
                                                                  $new_db2range->getPos2, 
                                                                  $new_db2range->getStrand);

  my $pseudo_sliceDB1 = $self->{'sliceadaptor1'}->fetch_by_region($new_db1range->getSeqLevel,
                                                                  $new_db1range->getSeqRegion,
                                                                  $new_db1range->getPos1, 
                                                                  $new_db1range->getPos2,
                                                                  $new_db1range->getStrand);

  my $pseudoDB2LoGene = new Bio::EnsEMBL::Gene(-DESCRIPTION => 'pseudo',
                                               -STABLE_ID   => 'pseudo',
                                               -SLICE       => $pseudo_sliceDB2,
                                               -START       => 1,
                                               -END         => $pseudo_sliceDB2->length,
                                               -STRAND      => 1);

  my $pseudoDB1LoGene = new Bio::EnsEMBL::Gene(-DESCRIPTION => 'pseudo',
                                               -STABLE_ID   => 'pseudo',
                                               -SLICE       => $pseudo_sliceDB1,
                                               -START       => 1,
                                               -END         => $pseudo_sliceDB1->length,
                                               -STRAND      => 1);

  my $new_boundary = new protein_boundary_ens(-ORI_SLICE   => $oldBoundary->ori_slice,
                                              -GENE1_1     => $oldBoundary->DB1HiGene,
                                              -GENE1_2     => $pseudoDB1LoGene,
                                              -GENE2_1     => $oldBoundary->DB2HiGene,
                                              -GENE2_2     => $pseudoDB2LoGene,
                                              -COMPLICATED => 1,
                                              -DB1DESC        => "ExtendDB1Forward"
                                             );
 
  return $new_boundary;
  #print "The direction is $guessed_direction_for_extension\n";
  #print "The distance guessed is $guessed_distance\n";
  #print "The new genomedb is ".$new_dbrange->getSequenceName."\n";
}

sub guessForwardDB2Direction {
  my $loweridx = shift;
  my $ordered_raw_genome_db_ranges_ref = shift;
  my $homolog_pairs_ref1 = shift;
  my $lower_raw_genomedb1range = $ordered_raw_genome_db_ranges_ref->[$loweridx];
  my $lower_genomedb2range     = new genomedbrange(-dbrange=>$homolog_pairs_ref1->{$lower_raw_genomedb1range}->[1]->feature_Slice->name);
  
  my $forwardCount  = 0;
  my $backwardCount = 0;

  for (my $i = 1; $i < 4; $i++) {
     my $scanIdx = $loweridx - $i;
     if ($scanIdx < 0) {last;}
     my $current_raw_genomedb1range = $ordered_raw_genome_db_ranges_ref->[$scanIdx];
     my $current_genomedb2range     = new genomedbrange(-dbrange=>$homolog_pairs_ref1->{$current_raw_genomedb1range}->[1]->feature_Slice->name);

     if ($current_genomedb2range->getSeqRegion ne $lower_genomedb2range->getSeqRegion) {next;}
     if ($current_genomedb2range->compare($lower_genomedb2range) > 0) {
        $backwardCount++;
     } else {
        $forwardCount++;
     }
  }

  if ($backwardCount > $forwardCount) {
     return -1;
  } elsif ($backwardCount < $forwardCount) {
     return 1;
  } else {
     warn "very difficult in forward near ".$lower_raw_genomedb1range."\n";
     return 0;
  }
}

sub guessDistance {
  my $genomedb1 = shift;
  my $genomedb2 = shift;
  my $factor = shift;
  if (!defined($factor)) {$factor = 3;}

  return ($genomedb1->getDistance($genomedb2) * $factor);
}

sub guessBackwardDB2Direction {
  my $hi_idx = shift;
  my $ordered_raw_genome_db_ranges_ref = shift;
  my $homolog_pairs_ref1 = shift;
  my $hi_raw_genomedb1range = $ordered_raw_genome_db_ranges_ref->[$hi_idx];
  my $hi_genomedb2range     = new genomedbrange(-dbrange=>$homolog_pairs_ref1->{$hi_raw_genomedb1range}->[1]->feature_Slice->name);
  
  my $forwardCount  = 0;
  my $backwardCount = 0;
  for (my $i = 1; $i < 4; $i++) {
     my $scanIdx = $hi_idx + $i;
     if ($scanIdx > scalar(@$ordered_raw_genome_db_ranges_ref)) {last;}
     my $current_raw_genomedb1range = $ordered_raw_genome_db_ranges_ref->[$scanIdx];
     my $current_genomedb2range     = new genomedbrange(-dbrange=>$homolog_pairs_ref1->{$current_raw_genomedb1range}->[1]->feature_Slice->name);

     if ($current_genomedb2range->getSeqRegion ne $hi_genomedb2range->getSeqRegion) {next;}
     if ($current_genomedb2range->compare($hi_genomedb2range) > 0) {
        $forwardCount++;
     } else {
        $backwardCount++;
     }
  }

  if ($backwardCount > $forwardCount) {
     return -1;
  } elsif ($backwardCount < $forwardCount) {
     return 1;
  } else {
     warn "very difficult in backward near ".$hi_raw_genomedb1range."\n";
     return 0;
  }
}

sub extendDB2Boundary {
  my $self = shift;
  my $oldBoundary = shift;
  my $target_genomedbrange = shift;
  
  my $homolog_pairs_ref1 = $self->{'homolog_pairs_ref1'};

  my $hi_genomedb1range = new genomedbrange(-dbrange=>$oldBoundary->DB1HiGene->feature_Slice->name);

  my $hi_genomedb2range = new genomedbrange(
     -dbrange=>$homolog_pairs_ref1->{$hi_genomedb1range->getSequenceName}->[1]->feature_Slice->name);

  my $low_genomedb1range = new genomedbrange(-dbrange=>$oldBoundary->DB1LoGene->feature_Slice->name);

  my $low_genomedb2range = new genomedbrange(
     -dbrange=>$homolog_pairs_ref1->{$low_genomedb1range->getSequenceName}->[1]->feature_Slice->name);

  my $ordered_raw_genome_db_ranges_by_db1_ref = $self->{'ordered_raw_genome_db_ranges_by_db1'};

  my $hi_idx  = dividefinder(0,                  scalar(@$ordered_raw_genome_db_ranges_by_db1_ref)-1, 
                             $hi_genomedb1range, $ordered_raw_genome_db_ranges_by_db1_ref);

  my $low_idx = dividefinder(0,                   scalar(@$ordered_raw_genome_db_ranges_by_db1_ref)-1,
                             $low_genomedb1range, $ordered_raw_genome_db_ranges_by_db1_ref);
  
  my $guessed_forward_direction_for_extension = 
      guessForwardDB2Direction($hi_idx, 
                               $ordered_raw_genome_db_ranges_by_db1_ref,
                               $homolog_pairs_ref1
                              );

  my $guessed_backward_direction_for_extension = -1 * 
      guessBackwardDB2Direction($low_idx,
                                $ordered_raw_genome_db_ranges_by_db1_ref,
                                $homolog_pairs_ref1
                              );

  my $guessed_forward_distance  = guessDistance($hi_genomedb1range,  $target_genomedbrange);
  my $guessed_backward_distance = guessDistance($low_genomedb1range, $target_genomedbrange);

  
  my $new_forward_db2range      = $hi_genomedb2range->shiftDistance($guessed_forward_distance, 
                                                                    $guessed_forward_direction_for_extension);
  
  #print $guessed_backward_direction_for_extension, "\n";
  #print $guessed_backward_distance, "\n";

  my $new_backward_db2range     = $low_genomedb2range->shiftDistance($guessed_backward_distance,
                                                                     $guessed_backward_direction_for_extension);

  
  my $pseudo_forward_sliceDB2 = $self->{'sliceadaptor2'}->fetch_by_region(
                                                           $new_forward_db2range->getSeqLevel, 
                                                           $new_forward_db2range->getSeqRegion, 
                                                           $new_forward_db2range->getPos1, 
                                                           $new_forward_db2range->getPos2, 
                                                           $new_forward_db2range->getStrand);

  my $pseudo_backward_sliceDB2 = $self->{'sliceadaptor2'}->fetch_by_region(
                                                           $new_backward_db2range->getSeqLevel, 
                                                           $new_backward_db2range->getSeqRegion, 
                                                           $new_backward_db2range->getPos1, 
                                                           $new_backward_db2range->getPos1 + 1, 
                                                           $new_backward_db2range->getStrand);

  my $pseudo_forward_DB2Gene = new Bio::EnsEMBL::Gene(-DESCRIPTION => 'pseudo',
                                                      -STABLE_ID   => 'pseudo',
                                                      -SLICE       => $pseudo_forward_sliceDB2,
                                                      -START       => 1,
                                                      -END         => $pseudo_forward_sliceDB2->length,
                                                      -STRAND      => 1);
  
  my $pseudo_backward_DB2Gene = new Bio::EnsEMBL::Gene(-DESCRIPTION => 'pseudo',
                                                       -STABLE_ID   => 'pseudo',
                                                       -SLICE       => $pseudo_backward_sliceDB2,
                                                       -START       => 1,
                                                       -END         => $pseudo_backward_sliceDB2->length,
                                                       -STRAND      => 1);

  
  my $new_forward_boundary = new protein_boundary_ens(-ORI_SLICE   => $oldBoundary->ori_slice,
                                                      -GENE1_1     => $oldBoundary->DB1HiGene,
                                                      -GENE1_2     => $oldBoundary->DB1LoGene,
                                                      -GENE2_1     => $oldBoundary->DB2HiGene,
                                                      -GENE2_2     => $pseudo_forward_DB2Gene,
                                                      -COMPLICATED => 1,
                                                      -DB1DESC     => $oldBoundary->getDB1Desc,
                                                      -DB2DESC     => "ExtendDB2Forward"
                                                     );

  
  my $new_backward_boundary = new protein_boundary_ens(-ORI_SLICE   => $oldBoundary->ori_slice,
                                                       -GENE1_1     => $oldBoundary->DB1HiGene,
                                                       -GENE1_2     => $oldBoundary->DB1LoGene,
                                                       -GENE2_1     => $pseudo_backward_DB2Gene,
                                                       -GENE2_2     => $oldBoundary->DB2LoGene,
                                                       -COMPLICATED => 1,
                                                       -DB1DESC     => $oldBoundary->getDB1Desc,
                                                       -DB2DESC     => "ExtendDB2Backward"
                                                      );

  return $new_forward_boundary, $new_backward_boundary;
}

sub getSliceAdp1 {
  my $self = shift;
  return $self->{'sliceadaptor1'};
}

sub getSliceAdp2 {
  my $self = shift;
  return $self->{'sliceadaptor2'};
}

1;
