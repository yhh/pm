package yenhua_arena;
use strict;
use Exporter();
use vars qw(
    @ISA
    @EXPORT
    @EXPORT_OK
    %EXPORT_TAGS
    $username
    $hostname
    $dbname
);
@ISA=qw(Exporter); 
@EXPORT=qw(
           $username $hostname $dbname
          ); 
$username='ensro'; 
$hostname='ecs1c'; 
$dbname='yenhua_arena'; 
1; 
