package reverse_transcribe;
use strict;
use Exporter();

our(
    @EXPORT,
    @ISA
);

@ISA=qw(Exporter);
@EXPORT=qw(&reverse_transcribe);

sub reverse_transcribe {
  my $seq = shift;
  my $newseq = $seq;
  if (ref($seq)) {
     $newseq = $$seq;
  }
  $newseq=~tr/AUTGCautgc/TAACGtaacg/;
  $newseq=reverse $newseq;
  return \$newseq;
}

1;
