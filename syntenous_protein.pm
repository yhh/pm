package syntenous_protein;
use strict;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::Root::Root;
use gene_rec;
use protein_boundary;

use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root);
$ID = 'syntenous_protein';
$VERSION  = 0.1;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($conf_file, $host, $dbuser, $port, $dbname, $assembly1, $assembly2, $species1, $species2) = 
  $self->_rearrange([qw(CONF_FILE HOST DBUSER PORT DBNAME ASSEMBLY1 ASSEMBLY2 SPECIES1 SPECIES2)], @args);
  
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
     $ensdb[$i] = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dbconf[$i]->{'host'},
                                                     -user   => $dbconf[$i]->{'user'},
                                                     -dbname => $dbconf[$i]->{'dbname'}
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

  $self->{'sliceadaptor'} = $db->get_db_adaptor($species1,$assembly1)->get_SliceAdaptor;

  my $gdadp        = $db->get_GenomeDBAdaptor;
  my $gdb1         = $gdadp->fetch_by_name_assembly($species1, $assembly1);
  my $gdb2         = $gdadp->fetch_by_name_assembly($species2, $assembly2);

  $self->{'genome_db_id1'} = $gdb1->dbID;
  $self->{'genome_db_id2'} = $gdb2->dbID;
  $self->{'homoloadp'}     = $db->get_HomologyAdaptor;

  my @seq_regs= @{$self->{'sliceadaptor'}->fetch_all('toplevel')};

  my %chr=();
  foreach my $seq_r(@seq_regs) {
    $chr{$seq_r->seq_region_name} = $seq_r;
  }

  $self->{'chr_regs'} = \%chr;

  return $self; # success - we hope!
}

sub getDB1SliceFromRange {
  my $self = shift;
  my $range= shift;
  if (!$range->isa('genomedbrange')) {
     warn "The range you provide is not a genomedbrange\n";
     return undef;
  }
  if ($range->getDBVersion ne $self->{'assembly1'}) {
     warn "The range you provide does not have the same db1 name\n";
     return undef;
  }
 return $self->ensDB1->get_SliceAdaptor->fetch_by_region($range->getSeqLevel, 
                                                         $range->getSequenceName,
                                                         $range->getPos1, 
                                                         $range->getPos2,
                                                         $range->getStrand
                                                        );
}

sub getDB2SliceFromRange {
  my $self = shift;
  my $range= shift;
  if (!$range->isa('genomedbrange')) {
     warn "The range you provide is not a genomedbrange\n";
     return undef;
  }
  if ($range->getDBVersion ne $self->{'assembly2'}) {
     warn "The range you provide does not have the same db1 name\n";
     return undef;
  }
 return $self->ensDB2->get_SliceAdaptor()->fetch_by_region($range->getSeqLevel, 
                                                           $range->getSequenceName,
                                                           $range->getPos1, 
                                                           $range->getPos2,
                                                           $range->getStrand
                                                          );
}

sub ensDB1 {
  my $self = shift;
  return $self->{'ensdb1'};
}

sub ensDB2 {
  my $self = shift;
  return $self->{'ensdb2'};
}

sub getBoundary {
  my $self     = shift;
  my $a_range  = shift;
  
  my $chr_regs= $self->{'chr_regs'};
  
  if (!$a_range->isa('simplerange')) {
     warn "$a_range is not a simplerange object\n";
     return undef;
  }
  
  if (!defined($chr_regs->{$a_range->getSeqName})) {
     warn $a_range->getSeqName." is not a valid chromosome\n";
     return undef;
  }

  my ($seq_region, $seq_region_start, $seq_region_end, $seq_region_strand) = 
     ($a_range->getSeqName,
      $a_range->getPos1,
      $a_range->getPos2, 
      $a_range->getStrand
     );

  my $slice = $self->{'sliceadaptor'}->fetch_by_region($chr_regs->{$seq_region}->coord_system->name,$seq_region);
  
  unless (defined $seq_region_start) {
    warn "WARNING : setting seq_region_start=1\n";
    $seq_region_start = 1;
  }

  if ($seq_region_start > $slice->length) {
     die "seq_region_start $seq_region_start larger than chr_length ".$slice->length."\n";
  }

  unless (defined $seq_region_end) {
     warn "WARNING : setting seq_region_end=seq_region->length ".$slice->length."\n";
     $seq_region_end = $slice->length;
  }

  if ($seq_region_end > $slice->length) {
     warn "WARNING : seq_region_end $seq_region_end larger than seq_region->length ".
     $slice->length."\n";
     $seq_region_end = $slice->length;
  }

  if (!defined($seq_region_strand)) {
     $seq_region_strand = 1;
  }

  if ($seq_region_start > $seq_region_end) {
     my $temp = $seq_region_start;
     $seq_region_start = $seq_region_end;
     $seq_region_end   = $temp;
     $seq_region_strand = -1;
  }

  my $ori_slice = $self->{'sliceadaptor'}->fetch_by_region('chromosome', $seq_region, $seq_region_start, $seq_region_end, $seq_region_strand);

  my $mem_ad    = $self->{'db'}->get_MemberAdaptor;
  $mem_ad->_final_clause("order by m.chr_start desc, m.chr_end desc limit 100");
  my @memberupperlist = sort by_chr_region_desc
     @{$mem_ad->_generic_fetch("m.chr_name = \'$seq_region\' and m.genome_db_id = ".$self->{'genome_db_id1'}.
                               " and m.chr_end < $seq_region_start")};

  $mem_ad->_final_clause("order by m.chr_start asc, m.chr_end asc limit 100");

  my @memberlowerlist = sort by_chr_region_asc 
      @{$mem_ad->_generic_fetch("m.chr_name = \'$seq_region\' and m.genome_db_id = ".$self->{'genome_db_id1'}.
                                " and m.chr_start > $seq_region_end")};
		
  my $upperGene =  $self->getNearestPCG(\@memberupperlist, 1);

  my $lowerGene =  $self->getNearestPCG(\@memberlowerlist, 1);

  my $upGene1 = $upperGene->[0]->{'-db1'};
  my $upGene2 = $upperGene->[0]->{'-db2'};
  my $loGene1 = $lowerGene->[0]->{'-db1'};
  my $loGene2 = $lowerGene->[0]->{'-db2'};

  if (defined($upGene1) && defined($upGene2) && defined($loGene1) && defined($loGene2)) {
     my $a_pair_of_boundary = new protein_boundary(-ORI_SLICE=>$ori_slice,
                                                   -GENE1_1=>$upGene1, 
                                                   -GENE2_1=>$upGene2,
                                                   -GENE1_2=>$loGene1,
                                                   -GENE2_2=>$loGene2,
                                                   -DB1    =>$self->ensDB1,
                                                   -DB2    =>$self->ensDB2
                                                  );
    return $a_pair_of_boundary; # should be objects of ensembl genes
  }
  return undef;
}


sub getNearestPCG {
	my $self = shift;
	my $refmemberlist  = shift;
	my $num_to_collect = shift;
	if (!defined($num_to_collect)) {
		$num_to_collect = 1;
	}
	my @memberlist = @$refmemberlist;
	#print "in getNearestPCG member:\t", scalar(@memberlist), "\n";
	my $count = 0;
	my @returnlist =();
	foreach my $single_mem (@memberlist) {
		my $description = $single_mem->description;
		my @homololist = @{$self->{'homoloadp'}->fetch_by_Member_paired_species($single_mem, $self->{'species2'})};
		foreach my $single_homolog (@homololist) {
			if ($single_homolog->description ne "UBRH") {next;}
			if ($count >= $num_to_collect) {last;}
			my @homologmemattr = @{$single_homolog->get_Member_Attribute_by_source("ENSEMBLGENE")};
			foreach my $RefArrayOfMemberAttributeArrayRef (@homologmemattr) {
				my $homolog_mem_attr = $RefArrayOfMemberAttributeArrayRef->[0];
				my $homolog_attr_attr= $RefArrayOfMemberAttributeArrayRef->[1];
				if ($single_mem->chr_name=~/NT/) {next;}
				if ($single_mem->stable_id eq $homolog_mem_attr->stable_id) {next;}
				# check chromosomes here?
				#
				$count++;
				my $gene1     = $self->{'genead1'}->fetch_by_stable_id($single_mem->stable_id);
				my $genename1 = $gene1->external_name;

				my $gene2     = $self->{'genead2'}->fetch_by_stable_id($homolog_mem_attr->stable_id);
				my $genename2 = $gene2->external_name;

				if (!defined($genename1)) {$genename1 = $gene1->display_id;}
				if (!defined($genename2)) {$genename2 = $gene2->display_id;}

                                #print join("\t", $genename1, $genename2), "\n";

				my $gene1= new gene_rec(-stable_id => $single_mem->stable_id,
							-gene_name => $genename1,
							-chr_name  => $single_mem->chr_name,
							-pos1      => $single_mem->chr_start,
							-pos2      => $single_mem->chr_end,
							-strand    => $single_mem->chr_strand,
							-source_name=>$single_mem->source_name
							);
				my $gene2= new gene_rec(-stable_id => $homolog_mem_attr->stable_id,
							-gene_name => $genename2,
							-chr_name  => $homolog_mem_attr->chr_name,
							-pos1      => $homolog_mem_attr->chr_start,
							-pos2      => $homolog_mem_attr->chr_end,
							-strand    => $homolog_mem_attr->chr_strand,
							-perc_id   => $homolog_attr_attr->perc_id,
							-perc_pos  => $homolog_attr_attr->perc_pos,
							-perc_cov  => $homolog_attr_attr->perc_cov
							 );
				my $dbgenes = {-db1=>$gene1, -db2=>$gene2};
				push(@returnlist, $dbgenes); 
			}
		}
	}
	return \@returnlist;
}


sub by_chr_region_desc {
    ($b->chr_start <=> $a->chr_start || $b->chr_end <=> $a->chr_end);
}

sub by_chr_region_asc {
    ($a->chr_start <=> $b->chr_start || $a->chr_end <=> $b->chr_end);
}

1;
