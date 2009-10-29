package test2;
use strict;

sub new {
  my $self = shift;
  bless {}, $self;
  return $self;
}

sub kk {
  my $self = shift;
  my ($test1, $test2) = (shift, shift);
  my @results = sort sortfunc2 ($test1, $test2);
  kk3("This is a test");
  return @results;
}

sub kk3 {
my $kk = shift;
print $kk,"\n";
}

sub sortfunc ($$) {
  print $_[0], "\n";
  print $_[1], "\n";
  return ($_[0] <=> $_[1]);
}

sub sortfunc2 {
  return ($a <=> $b);
}

1;
