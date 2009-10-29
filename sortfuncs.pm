package sortfuncs;
use strict;
use genomedbrange;
use string_hash_code;
use Exporter();

our(
    @EXPORT,
    @ISA
);

@ISA=qw(Exporter);
@EXPORT=qw(&orderbygenomedbrange &orderbymethoddbrange &orderbygenomedbrange_raw);

sub orderbygenomedbrange ($$) {
    my ($range1, $range2, @args) = @_;
    if ((!$range1->isa("genomedbrange")) || (!$range2->isa("genomedbrange"))) {
       die "The ranges you've provided are not valid genomedbranges, please check.\n";
    }

    return (hashcode($range1->getSeqRegion) <=> hashcode($range2->getSeqRegion) || 
            $range1->getPos1 <=> $range2->getPos1 || 
            $range1->getPos2 <=> $range2->getPos2);
}

sub orderbymethoddbrange ($$) {
    my ($obj1, $obj2, @args) = @_;

    my $range1 = $obj1->dbRange();
    my $range2 = $obj2->dbRange();

    if ((!$range1->isa("genomedbrange")) || (!$range2->isa("genomedbrange"))) {
       die "The ranges you've provided are not valid genomedbranges, please check.\n";
    }

    return (hashcode($range1->getSeqRegion) <=> hashcode($range2->getSeqRegion) || 
            $range1->getPos1 <=> $range2->getPos1 || 
            $range1->getPos2 <=> $range2->getPos2);
}

sub orderbygenomedbrange_raw ($$) {
    my ($range1_raw, $range2_raw, @args) = @_;
    my $range1 = new genomedbrange(-dbrange=>$range1_raw);
    my $range2 = new genomedbrange(-dbrange=>$range2_raw);
    
    if ((!$range1->isa("genomedbrange")) || (!$range2->isa("genomedbrange"))) {
       die "The ranges you've provided are not valid genomedbranges, please check.\n";
    }

    return (hashcode($range1->getSeqRegion) <=> hashcode($range2->getSeqRegion) || 
            $range1->getPos1 <=> $range2->getPos1 || 
            $range1->getPos2 <=> $range2->getPos2);
}
1;
