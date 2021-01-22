#!/usr/bin/perl -w
use strict;

my %fasta;
open FA, $ARGV[0] || die $!;
$/ = ">"; <FA>;
my $i = 1;
while(<FA>){
    chomp;
    $fasta{$i} = $_;
    $i++;
}
close FA;
$/ = "\n";

my %random = &GetRandNum($i, 5000);
foreach my $index (sort keys %random){
    my $out = $fasta{$index};
    print ">$out";
}

sub GetRandNum {
    my ($total, $num) = @_;

    my %hash;
    while ((keys %hash) < $num){
        $hash{int(rand{$total})} = 1;
    }

    return %hash;
}


__END__