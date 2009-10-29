package infernal;
use strict;
use Bio::Root::Root;
use Bio::Root::IO;
use vars qw($ID $VERSION @ISA);

@ISA= qw(Bio::Root::Root Bio::Root::IO);
$ID = 'infernal';
$VERSION  = 0.1;

=head1 

=cut

my $sequence_name=undef;
my $sequence_start=undef;
my $sequence_end=undef;

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  # initialize IO
  $self->_initialize_io(@args);

  return $self; # success - we hope!
}

sub next_sequence {
   my ($self) = @_;
   while (my $aline=$self->_readline()) {
     chomp($aline);
     #sequence: AC122206.4/171680-171904
     if ($aline=~/sequence:\s+(\S+)/) {
        $sequence_name=$1;
        if ($aline=~/sequence:\s+(\S*?)\/((\d+)\-(\d+))?/) {
           if ($2 ne '') {
              $sequence_start=$3;
              $sequence_end  =$4;
           } 
        }
        return $sequence_name;
     }
   }
   return undef;
}

sub next_hit {
  my ($self) = @_;
  my %onehit=();
  $onehit{'sequence_name'}=$sequence_name;
  $onehit{'sequence_start'}=$sequence_start;
  $onehit{'sequence_end'}=$sequence_end;

  {
    my $aline=$self->_readline();
    chomp($aline);
    if ($aline=~/^hit\s+(\d+)\s+:\s+(\d+)\s+(\d+)\s+(\S+)\s+bits/) {
       $onehit{'id'}=$1;
       $onehit{'hit_start'}=$2;
       $onehit{'hit_end'}=$3;
       $onehit{'bits'}=$4;
       last;
    } else {
       $self->_pushback($aline."\n");
       return undef;
    }
  } #after finding the word 'hit' as the starting point of reading infernal results
  
  while (my $aline=$self->_readline()) {
    chomp($aline);
    if ($aline=~/^sequence/ || $aline=~/^hit/ || $aline=~/CPU/) {
       $self->_pushback($aline."\n");
       last;
    }
    if ($aline eq '') {next;}

    # start read alignment here
    my $struc        =  $aline;
    my $std_seq      =  $self->_readline(); chomp($std_seq);
    my $CMstring     =  $self->_readline(); chomp($CMstring);
    my $target_seq   =  $self->_readline(); chomp($target_seq);

    if ($struc=~/\s+(\S+)\s*/) {
       $onehit{'struc'}.=$1;
    }
    
    my $prefix_space;
    if ($std_seq=~/(\s+)([\d|\-]+)(\s+)(\S+)\s*([\d|\-]+)/) {
       if (!defined($onehit{'std_start'})) {
          $onehit{'std_start'}=$2;
       }
       
       $onehit{'std_end'} =$5;
       $onehit{'std_seq'}.=$4;
    }
    $prefix_space = length($1.$2.$3);
    if ($CMstring=~/\s{$prefix_space}(.+)/) {
       $onehit{'CMstring'}.=$1;
    }
    
    if ($target_seq=~/(\s+)(\d+)(\s+)(\S+)\s*(\d+)/) {
       if (!defined($onehit{'target_start'})) {
          $onehit{'target_start'}=$2;
       }
       
       $onehit{'target_end'} =$5;
       $onehit{'target_seq'}.=$4;
    }
  }

  #process alignment here
  my @struc      = split(//, $onehit{'struc'});
  my @std_seq    = split(//, $onehit{'std_seq'});
  my @CMstring   = split(//, $onehit{'CMstring'});
  my @target_seq = split(//, $onehit{'target_seq'});

  my ($stem_count,     $cm_count, 
      $stem_offender,  $stem_conserved,    $stem_ambiguous,  
      $insertion,      $deletion)   = (0, 0, 0, 0, 0, 0, 0);
  
  for (my $i=0; $i<scalar(@struc); $i++) {
      if ($struc[$i] =~/[\(|\<|\>|\)]/) {
         $stem_count++;
         if    ($CMstring[$i] eq ' ') {$stem_offender++;}
         elsif ($CMstring[$i] eq ':') {$cm_count++;}
         elsif ($CMstring[$i] eq '+') {$stem_ambiguous++;}
         elsif ($CMstring[$i]=~/[A|U|G|C]/) {$stem_conserved++;}
         else  {die "Strange symbol in cm models: __$CMstring[$i]__ $sequence_name\n";}
      } elsif ($struc[$i] eq '.') {
         $insertion++;
      } 

      if ($target_seq[$i] eq '-') {
         $deletion++;
         
      }
  }
  
  $onehit{'stem_count'}     = $stem_count;
  $onehit{'cm_count'}       = $cm_count;
  $onehit{'stem_offender'}  = $stem_offender;
  $onehit{'stem_conserved'} = $stem_conserved;
  $onehit{'stem_ambiguous'} = $stem_ambiguous;
  $onehit{'insertion'}      = $insertion; # consider to count stem insertion and single-stranded insertion parts.
  $onehit{'deletion'}       = $deletion;  # consider to count the deletion part for stem and single-stranded parts.
  
  return \%onehit; 
}

1;
