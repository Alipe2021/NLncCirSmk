#!/usr/bin/perl -w
use strict;

###################################################################
#       Name: GetTPMFromSringtieGtfList.pl  
#       Auther: Liupeng
#       Email: sxliulian2012@hotmail.com
#       Date: 2021-01-18 08:17
####################################################################

die "Usage: perl $0 gtf_list.txt > out_tpm.matrix.tsv\n" if @ARGV < 1;

my %samples;
my %tpm_matrix;

open IN, $ARGV[0] || die $!;
while(<IN>){
    $_ =~ s/[\r\n]$//;    
    my ($sample, $gtf) = split /\t/, $_;
    $samples{$sample} = 1;
    
    open GTF, "<$gtf" || die "Can't open gtf file: $gtf\n";
    while(<GTF>){
        chomp;
        next if /^#/;
        my @a = split /\t/, $_;
        if ($a[2] eq "transcript"){
            my ($t_id) = $a[8] =~ /transcript_id "(\S+)";/;
            my ($tpm) = $a[8] =~ /TPM "(\S+)";/;
            $tpm_matrix{$t_id}{$sample} = $tpm;
        }else{
            next;
        }
    }
    close GTF;
}
close IN;


print join("\t", "Transcript", sort keys %samples)."\n";
foreach my $transcript (sort keys %tpm_matrix){
    my @out;
    foreach my $sam (sort keys %samples){
        my $tpm = exists $tpm_matrix{$transcript}{$sam} ? $tpm_matrix{$transcript}{$sam} : 0.0;
        push @out, $tpm;
    }

    print join("\t", $transcript, @out)."\n";
}

__END__
