package rfam;
use strict;
use Bio::Root::Root;
use DBI;

use vars qw($ID $VERSION @ISA);
@ISA= qw(Bio::Root::Root);
$ID = 'rfam';
$VERSION  = 0.1;


sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($host, $dbname, $dbuser, $pass, $port) = 
     $self->_rearrange([qw(HOST DBNAME DBUSER PASS PORT)], @args);

  if (!defined($port)) {$port = 3306;}
  if (!defined($host)) {$host = "";}
  if (!defined($dbuser)) {$dbuser = "yhh";}
  if (!defined($dbname)) {$dbname = "rfam8_1";}
  
  my $dbh = DBI->connect("DBI:mysql:dbname=".$dbname.";host=".$host.";port=".$port, $dbuser, $pass) || 
            die "#Cannot connect to DB $dbname using username $dbuser!\n";
  
  my $sth = $dbh->prepare('select t1.accession, t1.class, t1.id, t1.GA, t1.cmsearch from rfam_id t1');

  $sth->execute();

  my $rfam_id_hash_ref = {};
  my $rfam_acc_hash_ref = {};
  my $rfam_class_hash_ref = {};

  while (my $ref1 = $sth->fetchrow_hashref()) {
    my $oneRec = {
                  class    =>$ref1->{'class'},
                  id       =>$ref1->{'id'},
                  TC       =>$ref1->{'GA'},
                  acc      =>$ref1->{'accession'},
                  cmsearch =>$ref1->{'cmsearch'}
                 };
    
    $rfam_id_hash_ref->{$oneRec->{'id'}} = $oneRec;
    $rfam_acc_hash_ref->{$oneRec->{'acc'}} = $oneRec;

    my $recs_ref = $rfam_class_hash_ref->{$oneRec->{'class'}};
    if (!defined($recs_ref)) {$recs_ref = []; $rfam_class_hash_ref->{$oneRec->{'class'}} = $recs_ref;} 
    push(@$recs_ref, $oneRec);
  }

  $self->{'rfam_by_id'}    = $rfam_id_hash_ref;
  $self->{'rfam_by_acc'}   = $rfam_acc_hash_ref;
  $self->{'rfam_by_class'} = $rfam_class_hash_ref;

  return $self;
}

sub getTCByRfamName {
  my $self = shift;
  my $name = shift;
  return $self->{'rfam_by_id'}->{$name}->{'TC'};
}

sub getClassByRfamName {
  my $self = shift;
  my $name   = shift;
  return $self->{'rfam_by_id'}->{$name}->{'class'};
}

sub getCMSearchByRfamName {
  my $self     = shift;
  my $name     = shift;
  my $cmstring = $self->{'rfam_by_id'}->{$name}->{'cmsearch'};
  return $cmstring;
}

sub getAccByRfamName {
  my $self = shift;
  my $name = shift;
  return $self->{'rfam_by_id'}->{$name}->{'acc'};
}

sub getRfamNameByAcc {
  my $self = shift;
  my $acc  = shift;
  return $self->{'rfam_by_acc'}->{$acc}->{'id'};
}

sub getClassByAcc {
  my $self = shift;
  my $acc  = shift;
  return $self->{'rfam_by_acc'}->{$acc}->{'class'};
}

sub getCMSearchByAcc {
  my $self = shift;
  my $acc  = shift;
  return $self->{'rfam_by_acc'}->{$acc}->{'cmsearch'};
}

sub getTCByAcc {
  my $self = shift;
  my $acc  = shift;
  return $self->{'rfam_by_acc'}->{$acc}->{'TC'};
}

1;
