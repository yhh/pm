package math;
use strict;
use Exporter();

our(
    @EXPORT,
    @ISA
);

@ISA=qw(Exporter);
@EXPORT=qw(&min &max &mean &median &stdev &shuffle &shuffleCycles);

sub min {
    my $result;
    foreach my $single_value (@_) {
       if (defined($result)) {
          if ($result > $single_value) {$result = $single_value;}
       } else {
          $result = $single_value;
       }
    }
    return $result;
}

sub max {
    my $result;
    foreach my $single_value (@_) {
       if (defined($result)) {
          if ($result < $single_value) {$result = $single_value;}
       } else {
          $result = $single_value;
       }
    }
    return $result;
}

sub mean {
    my $result;
    foreach my $single_value (@_) {
       $result+=$single_value;
    }

    return $result/scalar(@_);
}

sub median {
    my $result; 
    my @list = sort @_;

    my $count = scalar(@list);
    
    if ($count % 2 == 0) {
       my $mid = $count / 2;
       $result = ($list[$mid] + $list[$mid +1]) / 2;
    } else {
       my $mid = int($count / 2);
       $result = $list[$mid];
    }
    return $result;
}

sub stdev {
    my @list = sort @_;

    my $mean = mean(@_);
    my $variance = 0;
    
    foreach my $single_value (@list) {
       my $difference = abs($single_value - $mean);
       $variance += $difference * $difference;
    }
    
    my $stderr = $variance / scalar(@list);
    my $stdev  = sqrt($stderr);

    return $stdev;
}

sub shuffle {
   my $list_ref = shift;
   return shuffleCycles($list_ref, 1000);
}

sub shuffleCycles {
  my $list_ref = shift;
  my $cycles   = shift;
 
  my @list = @$list_ref;
  
  for (my $i = 0; $i < $cycles; $i++) {
     my $idx1 = int(rand(@list));
     my $idx2 = int(rand(@list));
     my $temp = $list[$idx1];

     $list[$idx1] = $list[$idx2];
     $list[$idx2] = $temp;
  }

  return \@list;
}

1;
