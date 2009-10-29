package chr_format;
use strict;
use Bio::Root::Root;

our(@ISA);

@ISA=qw(Bio::Root::Root);

sub new {
  my $self = shift;
  my $region = shift;

  if ($region=~/(.*?)\/(\d+)\-(\d+)(\.(.))?/) {
     my $chr_string = $1;
     my $pos1       = $2;
     my $pos2       = $3;
     my $strand_str = $5;
     my $chr_name   = undef;

     if ($chr_string=~/chr(.*)/) {
        $chr_name   = $1; 
     } else {
        $chr_name   = $chr_string;
     }

     my $strand = undef;
     if (!defined($strand_str)) {
        $strand = 1;
     } elsif ($strand_str eq '-') {
        $strand = -1;
     } else {
        $strand = 1;
     }

     return bless {'chr_name' => $chr_name, 'start' => $pos1, 'end' => $pos2, 'strand' => $strand}, $self;

  } else {
     return undef;
  }
}

sub newChrPos {

  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  #$self->_initialize(@args);

  ($self->{'chr_name'}, $self->{'start'}, $self->{'end'}, $self->{'strand'}) = 
      $self->_rearrange([qw(CHR_NAME POS1 POS2 STRAND)], @args);
  
  if ($self->start > $self->end) {
     $self->{'strand'} = -1;
     my $temp = $self->{'start'};
     $self->{'start'} = $self->end;
     $self->{'end'}   = $temp;
  }
  
  if (!defined($self->{'strand'})) {
     $self->{'strand'} = 1;
  }
  return $self;
}

sub getChrString {
  my $self = shift;
  my $aLine = 'chr'.$self->chr_name.'/'.$self->start.'-'.$self->end;
  if ($self->strand == -1) {
     $aLine = $aLine.'.-';
  }
  return $aLine;
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

sub strand {
  my $self = shift;
  return $self->{'strand'};
}

1;
