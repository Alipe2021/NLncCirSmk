#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use vars qw($opt_b $opt_j $opt_o);

getopts('b:j:');
die "Usage: perl $0 -b bsj_matrix.tsv -j junc_ratio_matrix.tsv\n" if (! defined $opt_b || ! defined $opt_j);

my %hash;
open BSJ, $opt_b || die $!;
while(<BSJ>){
    chomp;
    next if /^$|^#/;
    next if /^CircRNA_ID/;

    my $n = 0;
    my @a = split /\t/, $_;
    for (my $i=1; $i<@a; $i++){
        $n++ if $a[$i] > 0;
    }
    if ($n > 1){
        $hash{$a[0]} = 1;
    }
}
close BSJ;

open JR, $opt_j || die $!;
while(<JR>){
    chomp;
    next if /^$|^#/;
    next if /^CircRNA_ID/;

    my @b = split /\t/, $_;
    my $m = 0;

    for (my $i=1; $i<@b; $i++){
        $m++ if $b[$i] >= 0.05;
    }
    if ($m > 0 && exists $hash{$b[0]}){
        print "$b[0]\n";
    }
}
close JR;


__END__
