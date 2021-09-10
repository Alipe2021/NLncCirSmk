#!/usr/bin/perl -w
use strict;

open IN, $ARGV[0] || die $!;
my %list = map{$_ =~ s/[\r\n]//g; $_, 1;}<IN>;
close IN;

open IN, $ARGV[1] || die $!;
my $header = <IN>; print $header;

while(<IN>){
    chomp;
    next if /^#|^$/;
    my @a = split /\t/, $_;

    print "$_\n" if (exists $list{$a[0]});
}
close IN;


__END__
