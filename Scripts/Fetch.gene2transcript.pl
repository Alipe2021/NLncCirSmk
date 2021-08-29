#!/usr/bin/env perl -w
use strict;

my %hash;
open GTF, $ARGV[0] || die $!;
while(<GTF>){
    chomp;
    next if /^$|^#/;
    my @a = split /\t/, $_;
    next if $a[2] ne "transcript";

    my ($transcript_id, $gene_id) = $a[8] =~ /transcript_id "(\S+)"; gene_id "(\S+)";/;
    $hash{$gene_id}{$transcript_id} = 1;
}
close GTF;

foreach my $g (sort keys %hash){
    foreach my $t (sort keys %{$hash{$g}}){
        print join("\t", $g, $t)."\n";
    }
}

__END__

