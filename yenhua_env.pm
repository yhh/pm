package yenhua_env;
use strict;
use Exporter();
use vars qw(
    @ISA
    @EXPORT
    @EXPORT_OK
    %EXPORT_TAGS
    $PRG
    $DATA
    $HOME
);
@ISA=qw(Exporter); 
@EXPORT=qw(
           $PRG $DATA
          ); 
$HOME='/nfs/team71/phd/yhh';
$PRG=$HOME.'/prg';
$DATA=$HOME.'/data';
1; 
