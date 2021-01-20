#!/usr/bin/perl -w
use strict;

open IN, $ARGV[0] || die $!;
my %list = map {chomp; $_, 1;}<IN>;
close IN;

my @files = glob("DESeq2.*.DEGs-FC4.xls");

my %hash;
foreach my $file (sort @files){
    my ($compared_name) = $file =~ /DESeq2.(\S+-vs-\S+).DEGs-FC4.xls/;
    open IN, $file || die $!;
    while(<IN>){
        chomp;
        my @a = split /\t/, $_;
        next if /^GeneID/;
        if (exists $list{$a[0]}){
            if ($a[2] > 0){
                $hash{$compared_name}{"UP"}++;
            }else{
                $hash{$compared_name}{"DOWN"}++;
            }
        }else{
            next;
        }
    }
    close IN;
}

print "Control-vs-Treatment\tUP\tDOWN\n";
foreach my $com (sort keys %hash){
    print join("\t", $com, $hash{$com}{"UP"}, $hash{$com}{"DOWN"})."\n";
}


__END__
