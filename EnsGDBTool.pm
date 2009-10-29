package EnsGDBTool;
use strict;
use Bio::Root::Root;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use genomedbrange;
use Bio::EnsEMBL::Registry;

use vars qw($ID $VERSION @ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT =
  qw(&getGDBAdpByConf &getSliceAdpByConf &checkOverlapByRawGDBR &dividefinder &getSliceAdpByReg);

$ID      = 'EnsGDBTool';
$VERSION = 0.1;

sub getGDBAdpByConf {
	my $conf_file = shift;

	#  print $conf_file,"\n";
	my $server_conf = do $conf_file
	  || die "Cannot find $conf_file or config file err!\n";

	my ( $host, $dbname, $user, $port ) = (
							   $server_conf->{'host'}, $server_conf->{'dbname'},
							   $server_conf->{'user'}, $server_conf->{'port'}
	);

	if ( $host   eq '' ) { die "Cannot find host $host!\n"; }
	if ( $dbname eq '' ) { die "Cannot find db $dbname!\n" }
	if ( $user   eq '' ) { die "Cannot find user $user!\n"; }
	if ( !defined($port) ) { $port = 3306; }

	my $dbadp =
	  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
										   -user   => $user,
										   -dbname => $dbname,
										   -host   => $host,
										   -port   => $port
	  );

	return $dbadp;
}

sub dividefinder {
	my ( $lower, $upper, $genomedbrange, $genomedbrangesRef ) =
	  ( shift, shift, shift, shift );

	if ( $upper - $lower == 1 ) {
		if ( $genomedbrange->compareRaw( $genomedbrangesRef->[$upper] ) > 0 ) {
			return $upper + 1;
		}
		elsif ( $genomedbrange->compareRaw( $genomedbrangesRef->[$lower] ) < 0 )
		{
			return $lower - 1;
		}
		else {
			return $lower;
		}
	}

	my $middle = int( ( $upper + $lower ) / 2 );

	#print join("\t", $upper, $middle, $lower), " upper middle lower\n";
	if (   ( $genomedbrange->compareRaw( $genomedbrangesRef->[$lower] ) >= 0 )
		&& ( $genomedbrange->compareRaw( $genomedbrangesRef->[$middle] ) < 0 ) )
	{
		return dividefinder( $lower, $middle, $genomedbrange,
							 $genomedbrangesRef );
	}
	else {
		return dividefinder( $middle, $upper, $genomedbrange,
							 $genomedbrangesRef );
	}
}

sub getSliceAdpByConf {
	my $conf_file = shift;
	my $dbadp     = getGDBAdpByConf($conf_file);

	if ( defined($dbadp) ) {
		return $dbadp->get_SliceAdaptor();
	}
	else {
		return undef;
	}
}

sub getSliceAdpByReg {
	my $registry = 'Bio::EnsEMBL::Registry';

	$registry->load_registry_from_db( -host => 'ensembldb.ensembl.org',
									  -user => 'anonymous',
									  -db_version => 54);
	return $registry->get_adaptor( 'Human', 'Core', 'Slice' );
}

sub checkOverlapByRawGDBR {
	my $raw_gdbr1 = shift;
	my $raw_gdbr2 = shift;

	my $gdbr1 = new genomedbrange( -dbrange => $raw_gdbr1 );
	my $gdbr2 = new genomedbrange( -dbrange => $raw_gdbr2 );

	if ( !defined($gdbr1) || !defined($gdbr2) ) {
		die
"Something wrong with the input raw ranges: $raw_gdbr1 or $raw_gdbr2\n";
	}

	return $gdbr1->checkOverlap($gdbr2);
}

1;
