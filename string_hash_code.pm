package string_hash_code;

use strict;
use Exporter();

our(
    @EXPORT,
    @ISA
);

@ISA=qw(Exporter);
@EXPORT=qw(&hashcode);

sub hashcode {
  my $result = '';
  my $temp_str = shift;
  if ($temp_str!~/\D+/) {return $temp_str;}
  my @fields = split(//, $temp_str);
  foreach my $single_asc (@fields) {
     $result = $result.ord($single_asc);
  }
  $result = $result + 0 ;
  return $result;
}

1;
