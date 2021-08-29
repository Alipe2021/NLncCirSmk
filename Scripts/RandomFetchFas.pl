#!/usr/bin/perl -w
use strict;

die "Usage: perl $0 in.fa [numbers] \n" if @ARGV < 1;

my $max_num = `grep -c ">" $ARGV[0]`;
my $numbers = $ARGV[1] ||= 1000;
my %random_num = &random($numbers, $max_num);

open FA, $ARGV[0] || die  $!;
$/ = ">";<FA>;
while(<FA>){
    chomp;
    if (exists $random_num{$.} ){
        print ">$_";
    }else{
        next;
    }
}
close FA;


##
sub random {
    my ($num, $max) = @_;
    my %hash;
    while ((keys %hash) < $num) {
        $hash{int(rand($max))} = 1;
    }
    return %hash;
}

__END__