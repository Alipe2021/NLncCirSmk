#!/usr/bin/perl -w
use strict;

## fetch transcript gtf 
open IN, $ARGV[0] || die $!;
my %hash = map {chomp; $_, 1;}<IN>;
close IN;

open GTF, $ARGV[1] || die $!;
while(<GTF>){
    chomp;
    if (/^#|^$/){
        print "$_\n";
        next;
    }else{
        my @a = split /\t/, $_;
        my ($id) = $a[8] =~ /transcript_id "(\S+)";/;

        if (exists $hash{$id}){
            print "$_\n";
        }else{
            next;
        }
    }    
}
close GTF;

__END__
