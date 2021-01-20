#!/usr/bin/perl -w
use strict;

my %length;
open IN, $ARGV[0] || die $!;
while(<IN>){
    chomp;
    my @a = split /\t/, $_;
    next if /Transcript/;

    if ($a[1] <= 200){
        $length{"<200"}++;
    }elsif(200 < $a[1] && $a[1] <= 500){
        $length{"200-500"}++;
    }elsif(500 < $a[1] && $a[1] <= 1000){
        $length{"500-1000"}++;
    }elsif(1000 < $a[1] && $a[1] <= 1500){
        $length{"1000-1500"}++;
    }elsif(1500 < $a[1] && $a[1] <= 2000){
        $length{"1500-2000"}++;
    }elsif(2000 < $a[1] && $a[1] <= 2500){
        $length{"2000-2500"}++;
    }elsif(2500 < $a[1] && $a[1] <= 3000){
        $length{"2500-3000"}++;
    }elsif(3000 < $a[1] && $a[1] <= 3500){
        $length{"3000-3500"}++;
    }elsif(3500 < $a[1] && $a[1] <= 4000){
        $length{"3500-4000"}++;
    }elsif(4000 < $a[1] && $a[1] <= 4500){
        $length{"4000-4500"}++;
    }elsif(4500 < $a[1] && $a[1] <= 5000){
        $length{"4500-5000"}++;
    }else{
        $length{">5000"}++;
    }

    if ($a[2] > 10){
        $length{">10"}++;
    }else{
        $length{$a[2]}++;
    }
}
close IN;

foreach my $key (sort keys %length){
    print "$key\t$length{$key}\n";
}

__END__
Transcript      Length  ExonNumber
ZeamMp002       348     1
ZeamMp003       1335    1
ZeamMp004       1164    1
ZeamMp005       2205    1