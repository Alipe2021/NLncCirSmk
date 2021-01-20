#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($intersection, $exclusion, $output);
GetOptions (
    "i|intersection=s" => \$intersection, 
    "e|exclusion:s"    => \$exclusion, 
    "o|output:s"  => \$output);

=head1 USGAE
    perl $0 --intersection a.txt,b.txt,c.txt --exclusion e1.txt,e2.txt,e3.txt -o out.list \n"
=cut

my @files = split /,/, $intersection;

my %hash;
foreach my $file (@files){
    open IN, $file || die $!;
    while(<IN>){
        chomp;
        $hash{$_}++;
    }
    close IN;
}

if (defined $exclusion){
    my @excl_files = split /,/, $exclusion;
    foreach my $ff (@excl_files){
        open IN, $ff || die $!;
        while(<IN>){
            chomp;
            delete $hash{$_};
        }
        close IN;
    }
}

open OUT, ">$output" || die $!;
foreach my $key (sort keys %hash){
    print OUT "$key\n" if $hash{$key} > 1;
}
close OUT;

__END__

