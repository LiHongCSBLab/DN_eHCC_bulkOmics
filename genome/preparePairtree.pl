
=pod
ARGV[0]: patientID
ARGV[1]: samples IDs, seperated by "," (The 1st sample is the control sample, and other samples are the disease samples),

=cut

my $minimumDepth = 15;
$home="/WGS/";
my $dir = "$home/6-SubcloneAnalysis/pairtree/muTect_varScan_overlap_clustervarsByPairtree/test_$ARGV[0]_minimumDepth$minimumDepth";
if( !-e $dir ){
  system(" mkdir $dir ");
}
my @samples = split(",", $ARGV[1]);

my %filter = ();
my $snvDir="$home/5-SomaticMutation/compare_muTect_varScan/".$ARGV[0].".list.wgs.overlap.readN/";
for(my $i=0; $i<@samples; $i++){
  my $sample = $samples[$i];
  open(In, "$snvDir/$sample.snv.readN") or die "Can not open $snvDir/$sample.snv.readN\n";
  while(<In>){
    chomp;
    my @array = split(/\t/, $_);
    if( ($array[0] eq "chrM") || ($array[0] eq "chrX") || ($array[0] eq "chrY") ){ print STDERR "$samples[$i] $_\n";; next; }
    $varReadN{"$array[0]:$array[1]:$array[2]>$array[3]"}{$sample} = $array[6];
    $refReadN{"$array[0]:$array[1]:$array[2]>$array[3]"}{$sample} = $array[5] - $array[6];
    if( $array[5] < $minimumDepth ){
      $filter{"$array[0]:$array[1]:$array[2]>$array[3]"} = 0;
      print STDERR "$samples[$i] $_\n";
    }
    if( $i==0 & $array[6]>0  ){
      $filter{"$array[0]:$array[1]:$array[2]>$array[3]"} = 0;
      print STDERR "$samples[$i] $_\n";
    }
  }
  close(In);
}

my $s = 0;
open(Out, ">$dir/$ARGV[0].ssm");
print Out "id\tname\tvar_reads\ttotal_reads\tvar_read_prob\n";
foreach my $key (sort keys %varReadN){
  if( !exists $filter{$key} ){
    $s = $s + 1;
    print Out "s$s\t$key";
    my @array1; my @array2; my @array3;
    for(my $i=1; $i<@samples; $i++){
      push @array1, $varReadN{$key}{$samples[$i]};
      push @array2, $varReadN{$key}{$samples[$i]}+$refReadN{$key}{$samples[$i]};
      push @array3, 0.5;
    }
    print Out "\t".join(",", @array1)."\t".join(",", @array2)."\t".join(",", @array3)."\n";
  }
}
close(Out);

open(Out, ">$dir/$ARGV[0].tmp.params.json");
print Out "{\"samples\": [\"".join("\", \"", @samples[1..(scalar(@samples)-1)])."\"], \"garbage\": []}\n";
close(Out);

system("/program/anaconda3/bin/python /program/pairtree-master/bin/clustervars $dir/$ARGV[0].ssm $dir/$ARGV[0].tmp.params.json $dir/$ARGV[0].params.json\n");
system("/program/pairtree-master/bin/pairtree --params $dir/$ARGV[0].params.json $dir/$ARGV[0].ssm $dir/$ARGV[0].pairtree.npz\n");
system("/program/pairtree-master/bin/plottree --runid $ARGV[0] $dir/$ARGV[0].ssm $dir/$ARGV[0].params.json $dir/$ARGV[0].pairtree.npz $dir/$ARGV[0].pairtree.html\n");
