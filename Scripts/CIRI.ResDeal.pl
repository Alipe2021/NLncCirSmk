#!/usr/bin/env perl -w
use strict;
use Data::Dumper;

die "Useage: perl $0 <GTF list> [output dir] \n" if @ARGV < 1;

open IN, $ARGV[0] || die "Can't open GTF list file: $$ARGV[0] \n";
my %fileHash = map{chomp; my @a = split /\t/, $_; @a;}<IN>;
close IN;

# print Dumper %fileHash;
my %InfoHash;
my %BsjHash;
my %JuncHash;
my %CpmHash;
my %SrpbmHash;

foreach my $sample (sort keys %fileHash){
    my $file = $fileHash{$sample};
    open GTF, $file || die "Can't open gtf file: $file\n Please check gtf files\n";
    while(<GTF>){
        chomp;
        next if /^$|^#/;
        my @a = split /\t/, $_;
        my ($circ_id) = $a[8] =~ /circ_id "(\S+:\d+\|\d+)";/;
        my ($circ_type) = $a[8] =~ /circ_type "(\w+)";/;
        my ($bsj) = $a[8] =~ /bsj (\d+\.\d+);/;
        my ($junc_ratio) = $a[8] =~ /junc_ratio (\d+\.\d+);/;

        my $gene_id;
        if (/gene_id/){
            ($gene_id) = $a[8] =~ /gene_id "(\S+)";/;
        }else{
            $gene_id = "--";
        }
        my $circ_len = $a[4] - $a[3] + 1;
        my $circ_srpbm = ($a[5] * 1000) / $circ_len;

        $InfoHash{$circ_id} = join("\t", $circ_id, @a[0,3,4,6], $circ_type, $circ_len, $gene_id);
        $BsjHash{$circ_id}{$sample} = $bsj;
        $JuncHash{$circ_id}{$sample} = $junc_ratio;
        $CpmHash{$circ_id}{$sample} = $a[5];
        $SrpbmHash{$circ_id}{$sample} = sprintf("%.3f", $circ_srpbm);
    }
    close GTF;
}
###
my $output_dir = $ARGV[1] ||= "./";
if (! -d $output_dir){
    mkdir $output_dir;
}
open O1, ">$output_dir/CircRNAs.Info.tsv" || die $!;
print O1 "CircRNA_ID\tChrom\tStart\tEnd\tStrand\tCircRNA_type\tCircRNA_Length\tHostGeneID\n";
foreach my $circ (sort keys %InfoHash){
    print O1 $InfoHash{$circ}."\n";
}
close O1;
###
open O2, ">$output_dir/CircRNAs.BSJ_Matrix.tsv" || die $!;
print O2 join("\t", "CircRNA_ID", sort keys %fileHash)."\n";
foreach my $c1 (sort keys %BsjHash){
    my @tmp;
    foreach my $sam (sort keys %fileHash){
        my $bsj = exists $BsjHash{$c1}{$sam} ? $BsjHash{$c1}{$sam} : "0.0";
        push @tmp, $bsj;
    }
    print O2 join("\t", $c1, @tmp)."\n";
}
close O2;
###
open O3, ">$output_dir/CircRNAs.JuncRatio_Matrix.tsv" || die $!;
print O3 join("\t", "CircRNA_ID", sort keys %fileHash)."\n";
foreach my $circ (sort keys %JuncHash){
    my @tmp;
    foreach my $sam (sort keys %fileHash){
        my $jr = exists $JuncHash{$circ}{$sam} ? $JuncHash{$circ}{$sam} : "0.0";
        push @tmp, $jr;
    }
    print O3 join("\t", $circ, @tmp)."\n";
}
close O3;
###
open O4, ">$output_dir/CircRNAs.CPM_Matrix.tsv" || die $!;
print O4 join("\t", "CircRNA_ID", sort keys %fileHash)."\n";
foreach my $circ (sort keys %CpmHash){
    my @tmp;
    foreach my $sam (sort keys %fileHash){
        my $cpm = exists $CpmHash{$circ}{$sam} ? $CpmHash{$circ}{$sam} : "0.0";
        push @tmp, $cpm;
    }
    print O4 join("\t", $circ, @tmp)."\n";
}
close O4;
###
open O5, ">$output_dir/CircRNAs.SRPBM_Matrix.tsv" || die $!;
print O5 join("\t", "CircRNA_ID", sort keys %fileHash)."\n";
foreach my $circ (sort keys %SrpbmHash){
    my @tmp;
    foreach my $sam (sort keys %fileHash){
        my $srpbm = exists $SrpbmHash{$circ}{$sam} ? $SrpbmHash{$circ}{$sam} : "0.0";
        push @tmp, $srpbm;
    }
    print O5 join("\t", $circ, @tmp)."\n";
}
close O5;

__END__
