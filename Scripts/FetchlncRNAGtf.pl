#!/usr/bin/perl -w
use strict;

my %hash;
open GTF, $ARGV[0] || die $!;
while(<GTF>){
    chomp;
    next if (/^#|^$/);
    my @a = split /\t/, $_;
    next if ($a[2] ne "transcript");

    my ($id) = $a[8] =~ /transcript_id "(\S+)";/;
    $hash{$id} = 1;
}
close GTF;

open GTF, $ARGV[0] || die $!;
while(<GTF>){
    chomp;
    next if /^#|^$/;
    my @a = split /\t/, $_;
    
    my ($id) = $a[8] =~ /transcript_id "(\S+)";/;
    if (exists $hash{$id}){
        print "$_\n";
    }else{
        next;
    }
}
close GTF;


__END__
1       StringTie       transcript      1207441 1213412 .       -       .       transcript_id "MSTRG.82.5"; gene_id "MSTRG.82"; xloc "XLOC_003174"; class_code "u"; tss_id "TSS6093";
1       StringTie       exon    1207441 1208003 .       -       .       transcript_id "MSTRG.82.5"; gene_id "MSTRG.82"; exon_number "1";
1       StringTie       exon    1208103 1208235 .       -       .       transcript_id "MSTRG.82.5"; gene_id "MSTRG.82"; exon_number "2";
1       StringTie       exon    1208762 1208789 .       -       .       transcript_id "MSTRG.82.5"; gene_id "MSTRG.82"; exon_number "3";
1       StringTie       exon    1208877 1209075 .       -       .       transcript_id "MSTRG.82.5"; gene_id "MSTRG.82"; exon_number "4";
1       StringTie       exon    1209861 1211495 .       -       .       transcript_id "MSTRG.82.5"; gene_id "MSTRG.82"; exon_number "5";
1       StringTie       exon    1212265 1213412 .       -       .       transcript_id "MSTRG.82.5"; gene_id "MSTRG.82"; exon_number "6";