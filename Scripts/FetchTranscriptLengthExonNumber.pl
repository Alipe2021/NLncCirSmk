#!/usr/bin/perl -w
use strict;

my %hash;
open GTF, $ARGV[0] || die $!;
while(<GTF>){
    chomp;
    next if /^#|^$/;
    my @a = split /\t/, $_;
    my ($transcript) = $a[8] =~ /transcript_id "(\S+)";/;
    
    if ($a[2] eq 'transcript'){
        $hash{$transcript}{'length'} = $a[4] - $a[3] + 1;
    }elsif ($a[2] eq 'exon'){
        $hash{$transcript}{'exon'}++;
    }else{
        next;
    }
}
close GTF;

print "Transcript\tLength\tExonNumber\n";
foreach my $id (sort keys %hash){
    print join("\t", $id, $hash{$id}{'length'}, $hash{$id}{'exon'})."\n";
}


__END__
1       StringTie       transcript      108856  109422  .       +       .       transcript_id "MSTRG.5.1"; gene_id "MSTRG.5"; xloc "XLOC_000002"; class_code "u"; tss_id "TSS2";
1       StringTie       exon    108856  109016  .       +       .       transcript_id "MSTRG.5.1"; gene_id "MSTRG.5"; exon_number "1";
1       StringTie       exon    109354  109422  .       +       .       transcript_id "MSTRG.5.1"; gene_id "MSTRG.5"; exon_number "2";
