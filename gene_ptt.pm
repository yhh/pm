package gene_ptt;

use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use gene_rec;

use vars qw($ID $VERSION @ISA);

@ISA     = qw(Bio::Root::Root  Bio::Root::IO);
$ID      = 'gene_ptt';
$VERSION = 0.1;

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	# initialize IO
	$self->_initialize_io(@args);

	return $self;    # success - we hope!
}

sub get_all_genes {
	my $self = shift;
	my @res  = ();
	while ( my $aGene = $self->next_gene ) {
		push( @res, $aGene );
	}

	if ( scalar(@res) == 0 ) {
		die "Cannot find any genes!\n";
	}

	my @sorted_res =
	  sort { $a->getPos1 <=> $b->getPos2 || $a->getStrand <=> $b->getStrand }
	  @res;

	return \@sorted_res;
}

sub next_gene {
	my $self = shift;
	while ( my $aline = $self->_readline() ) {
		chomp($aline);

		#2801..3733
		#+
		#310
		#-
		#thrB
		#b0003
		#-
		#homoserine kinase

		my @fields = split( "\t", $aline );
		if ( scalar(@fields) < 8 ) {
			next;
		}

		my (
			 $pos1,      $pos2,    $strand, $length, $PID,
			 $gene_name, $synonym, $COG,    $desc
		);

		if ( $fields[0] =~ /^(\d+)\.\.(\d+)/ ) {
			$pos1 = $1;
			$pos2 = $2;

			if ( $fields[1] eq '+' ) {
				$strand = 1;
			}
			else {
				$strand = -1;
			}

			$length    = $fields[2];
			$PID       = $fields[3];
			$gene_name = $fields[4];
			$synonym   = $fields[5];
			$COG       = $fields[6];
			$desc      = $fields[7];
		}
		else {
			next;
		}
		return
		  new gene_rec(
						-POS1      => $pos1,
						-POS2      => $pos2,
						-GENE_NAME => $gene_name,
						-STABLE_ID => $synonym,
						-STRAND    => $strand
		  );
	}
}

1;
