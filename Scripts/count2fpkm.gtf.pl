use strict;
use warnings;
die "Usage: perl $0 count.txt gene.gtf\n" if(scalar @ARGV !=2);
my %len;
open IN,$ARGV[1] or die $!; # gtf
while(<IN>){
    next if($_=~/#/);
    chomp;
    my($chr,undef,$type,$start,$end,undef,undef,undef,$str)=split /\t/,$_;
    next if($type ne "gene");
    my $gid;
    #if($str=~/gene_id\s+\"*([0-9a-zA-Z:._-]+)/){
    if($str=~/gene_id "(\S+)";/){
        $gid=$1;
    $gid=~s/\w+://g;
}else{
    die $str;
    die "GTF malformated\n";
}
$len{$gid}=abs($end-$start)+1;
}
close IN;

my %libsize;
open IN,$ARGV[0] or die $!; # count txt
my $header=<IN>;
print "$header";
my @header=split /,/,$header;
while(<IN>){
    chomp;
    my @line=split /,/,$_;
    for(my $i=1;$i<@line;$i++){
        $libsize{$header[$i]}+=$line[$i];
    }
}
seek IN,0,0;
<IN>;
while(<IN>){
    chomp;
    my @line=split /\s+/,$_;
    my $length=$len{$line[0]};
    print "$line[0]";
    for(my $i=1;$i<@line;$i++){
        my $fpkm=(1000000000)*$line[$i]/($length*$libsize{$header[$i]});
        if($fpkm>=0.01 or $fpkm==0){
            $fpkm=sprintf("%.2f",$fpkm);
        }else{
            $fpkm=sprintf("%.2e",$fpkm);
        }
        print ",$fpkm";
    }
#print "\t$length\n";
    print "\n";
}

