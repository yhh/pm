package gene_rec;
use strict;
use Bio::Root::Root;
use simplerange;

our (@ISA);

@ISA = qw(Bio::Root::Root);

sub new {

	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	#$self->_initialize(@args);

	(
	   $self->{'chr_name'},    $self->{'start'},     $self->{'end'},
	   $self->{'strand'},      $self->{'stable_id'}, $self->{'gene_name'},
	   $self->{'source_name'}, $self->{'perc_id'},   $self->{'perc_pos'},
	   $self->{'perc_cov'}
	  )
	  = $self->_rearrange(
		[
		  qw(CHR_NAME POS1 POS2 STRAND STABLE_ID GENE_NAME SOURCE_NAME PERC_ID PERC_POS PERC_COV)
		],
		@args
	  );

	return $self;
}

sub toSimpleRange {
	my $self = shift;
	return
	  new simplerange(
					   -POS1   => $self->getPos1,
					   -POS2   => $self->getPos2,
					   -DESC   => $self->gene_name,
					   -STRAND => $self->getStrand
	  );
}

sub toString2 {
	my $self = shift;
	return
	  join( "\t",
			$self->stable_id, $self->gene_name, $self->start,
			$self->end,       $self->strand );
}

sub toString3 {
	my $self = shift;
	return
	  join( "\t",
			$self->gene_name, $self->start,
			$self->end,       $self->strand );
}

sub chr_name {
	my $self = shift;
	return $self->{'chr_name'};
}

sub start {
	my $self = shift;
	return $self->{'start'};
}

sub end {
	my $self = shift;
	return $self->{'end'};
}

sub getPos1 {
	my $self = shift;
	return $self->{'start'};
}

sub getPos2 {
	my $self = shift;
	return $self->{'end'};
}

sub strand {
	my $self = shift;
	return $self->{'strand'};
}

sub getStrand {
	my $self = shift;
	return $self->{'strand'};
}

sub stable_id {
	my $self = shift;
	return $self->{'stable_id'};
}

sub gene_name {
	my $self = shift;
	return $self->{'gene_name'};
}

sub source_name {
	my $self = shift;
	return $self->{'source_name'};
}

sub perc_id {
	my $self = shift;
	return $self->{'perc_id'};
}

sub perc_pos {
	my $self = shift;
	return $self->{'perc_pos'};
}

sub perc_cov {
	my $self = shift;
	return $self->{'perc_cov'};
}

1;
