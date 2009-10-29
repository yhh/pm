package genetotissue;
use strict;
use Bio::Root::Root;
use DBI;

use vars qw($ID $VERSION @ISA);
@ISA= qw(Bio::Root::Root);
$ID = 'genetotissue';
$VERSION  = 0.1;


sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($host, $dbname, $dbuser, $pass, $port) = 
     $self->_rearrange([qw(HOST DBNAME DBUSER PASS PORT)], @args);

  if (!defined($port)) {$port = 3306;}
  if (!defined($host)) {$host = "";}
  if (!defined($dbuser)) {$dbuser = "yhh";}
  if (!defined($dbname)) {$dbname = "library";}
  
  my $dbh = DBI->connect("DBI:mysql:dbname=".$dbname.";host=".$host.";port=".$port, $dbuser, $pass) || 
            die "#Cannot connect to DB $dbname using username $dbuser!\n";
  
  my $sth = $dbh->prepare('select t1.gn_name, t1.LIB, t1.library, t1.tissue, t1.ug_tissue, t1.histology, t1.organism, t3.DEVELOPMENTAL_STAGE, t3.LID from gene2tissue t1, ug_alias_acc_lid t2, Hs_lib_info t3 where t1.library = ? and t2.ACC=t1.library and t3.LID = t2.LID');

  $self->{'sth'} = $sth;
  return $self;
}

sub getRecByAcc {
  my $self = shift;
  my $name = shift;
  my @res  = ();
  my $sth = $self->{'sth'};

  $sth->bind_param(1, $name);
  $sth->execute();

  while (my $ref1 = $sth->fetchrow_hashref()) {
    my $oneRec = {
                  gn_name  =>$ref1->{'gn_name'},
                  LIB      =>$ref1->{'LIB'},
                  library  =>$ref1->{'library'},
                  tissue   =>$ref1->{'tissue'},
                  ug_tissue=>$ref1->{'ug_tissue'},
                  histology=>$ref1->{'histology'},
                  organism =>$ref1->{'organism'},
                  DEVELOPMENTAL_STAGE =>$ref1->{'DEVELOPMENTAL_STAGE'},
                  LID      =>$ref1->{'LID'}
                 };

    push(@res, $oneRec);
  }
  return \@res;
}

1;
