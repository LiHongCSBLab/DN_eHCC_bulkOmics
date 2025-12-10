# perl compare_somatic.pl FFPE2.list

use strict;
use List::Util qw/min/;
my $home = "/WGS/";
my %wgs_snv; my %wgs_info;
my %snv; my %info;

open(In, $ARGV[0]) or die "Can not open $ARGV[0]\n";
while(<In>){
  chomp;
  my ($name, $file) = split(/\t/, $_);
  open(In1, "$home/$file") or die "Can not open $home/$file\n";
  <In1>;
  while(<In1>){
    chomp;
    my @array = split(/\t/, $_);
    if( $array[3] eq $array[4] ){ next; }
    if( $name=~/WGS/ ){ 
      $wgs_info{$_}{$name} = 1; 
      $wgs_snv{$name}{"$array[0]\t$array[1]\t$array[3]\t$array[4]"} = 1;
    }
    if( ($array[5] eq "exonic") || ($array[5] eq "splicing") ){
      $info{$_}{$name} = 1;
      $snv{$name}{"$array[0]\t$array[1]\t$array[3]\t$array[4]"} = 1;
    }
  }
  close(In1);
}
close(In);

my %snvN;
open(Out0, ">$ARGV[0].exonic.stat");
print Out0 "sample\tsnv_number\n";
foreach my $name (sort keys %snv){
  print Out0 "$name\t".scalar(keys %{$snv{$name}})."\n";
  $snvN{$name} = scalar(keys %{$snv{$name}});
}
close(Out0);

open(Out0, ">$ARGV[0].exonic.compare");
print Out0 "sample1\tsample2\toverlap\tsample1_unique\tsample2_unique\tjaccardSimilarity\tsimilarity2\n";
my @names = sort keys %snv;
for(my $i=0; $i<@names-1; $i++){
  for(my $j=$i+1; $j<@names; $j++){
    my $overlap = 0;
    my @snvs1 = keys %{$snv{$names[$i]}};
    foreach my $snv( @snvs1 ){
      if( exists $snv{$names[$j]}{$snv} ){
        $overlap = $overlap + 1;
      }
    }
    my $unique1 = $snvN{$names[$i]} - $overlap;
    my $unique2 = $snvN{$names[$j]} - $overlap;
    # Jaccard similarity coefficient
    my $similarity = $overlap/($overlap+$unique1+$unique2);
    my $similarity2 = $overlap/min($overlap+$unique1, $overlap+$unique2);
    print Out0 "$names[$i]\t$names[$j]\t$overlap\t$unique1\t$unique2\t$similarity\t$similarity2\n";
  }
}
close(Out0);

open(Out0, ">$ARGV[0].exonic.union");
foreach my $key (sort keys %info){
  print Out0 join(",",sort keys %{$info{$key}})."\t$key\n";
}
close(Out0);

my %wgs_snvN;
open(Out0, ">$ARGV[0].wgs.stat");
print Out0 "sample\tsnv_number\n";
foreach my $name (sort keys %wgs_snv){
  print Out0 "$name\t".scalar(keys %{$wgs_snv{$name}})."\n";
  $wgs_snvN{$name} = scalar(keys %{$wgs_snv{$name}});
}
close(Out0);

open(Out0, ">$ARGV[0].wgs.compare");
print Out0 "sample1\tsample2\toverlap\tsample1_unique\tsample2_unique\tjaccardSimilarity\tsimilarity2\n";
my @names = sort keys %wgs_snv;
for(my $i=0; $i<@names-1; $i++){
  for(my $j=$i+1; $j<@names; $j++){
    my $overlap = 0;
    my @wgs_snvs1 = keys %{$wgs_snv{$names[$i]}};
    foreach my $wgs_snv( @wgs_snvs1 ){
      if( exists $wgs_snv{$names[$j]}{$wgs_snv} ){
        $overlap = $overlap + 1;
      }
    }
    my $unique1 = $wgs_snvN{$names[$i]} - $overlap;
    my $unique2 = $wgs_snvN{$names[$j]} - $overlap;
    # Jaccard similarity coefficient
    my $similarity = $overlap/($overlap+$unique1+$unique2);
    my $similarity2 = $overlap/min($overlap+$unique1, $overlap+$unique2);
    print Out0 "$names[$i]\t$names[$j]\t$overlap\t$unique1\t$unique2\t$similarity\t$similarity2\n";
  }
}
close(Out0);

open(Out0, ">$ARGV[0].wgs.union");
foreach my $key (sort keys %wgs_info){
  print Out0 join(",",sort keys %{$wgs_info{$key}})."\t$key\n";
}
close(Out0);

=pod
grep "MuTect.*VarScan" $ARGV[0].wgs.union > test
grep "VarScan.*Mutect" $ARGV[0].wgs.union >> test
# head -1 $ARGV[0].wgs.union > $ARGV[0].wgs.overlap
sort test | uniq >> $ARGV[0].wgs.overlap
grep -v -P "\tChr\t" $ARGV[0].wgs.overlap | cut -f2,3,5,6  > $ARGV[0].wgs.overlap.bed

grep -P "\texonic\t" $ARGV[0].wgs.overlap > $ARGV[0].exonic.overlap
grep -P "\tsplicing\t" $ARGV[0].wgs.overlap >> $ARGV[0].exonic.overlap

=cut

