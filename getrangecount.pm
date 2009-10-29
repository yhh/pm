package getrangecount;
use DBI;
use yenhua_as;

my ($tablename, $fieldname)=('', '');

init();
main();
final();

sub init {
if (scalar(@ARGV) != 2) {die "Please input tablename and field to count the range!\n";}
$tablename=shift @ARGV; $fieldname=shift @ARGV;

$dbh = DBI->connect("DBI:mysql:dbname=$dbname;host=$hostname", $username , undef) || 
       die "#Cannot connect to DB $dbname using username $username to retrieve!\n";
}


sub main {
    my $sth  =$dbh->prepare("select min($fieldname), max($fieldname) from $tablename") || 
               die "Cannot prepare SQL using field $fieldname from table $tablename!\n";
    my $sth2 =$dbh->prepare("select count(*) from $tablename where $fieldname >= ? and $fieldname <?") || 
               die "Cannot prepare SQL for intervals!\n";

    $sth->execute();
    my @range=$sth->fetchrow_array();
    if ($range[1] == 0) {die "The maximal value of the field $fieldname is zero. Please check it!\n";}
    
    my $start=0; my $end=1;
    while ($start < $range[1]) 
      {
       $sth2->bind_param(1, $start); $sth2->bind_param(2, $end);
       $sth2->execute(); my @count=$sth2->fetchrow_array(); 
       print "$start~$end\t$count[0]\n";
       $start=$end; $end=$end*2;
      }
    $sth->finish();
    $sth2->finish();
}

sub final {
    $dbh->disconnect;
}
