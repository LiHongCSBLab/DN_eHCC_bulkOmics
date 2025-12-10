
open(In, $ARGV[0]) or die;
open(Out0, ">$ARGV[0].annovar");
open(Out1, ">$ARGV[0].bed");
while(<In>){
  chomp;
  my @array = split(/\t/, $_);
  print Out0 "$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]\n";
  print Out1 "$array[1]\t$array[2]\t$array[4]\t$array[5]\n";
}
close(In);
close(Out0);
close(Out1);

# ANNOVAR
system("/meta/annovar/table_annovar.pl $ARGV[0].annovar /meta/annovar/humandb -buildver hg19 -protocol refGene,clinvar_20160302,cosmic70,snp138,esp6500siv2_all,ALL.sites.2015_08,kaviar_20150923,exac03nontcga,ljb26_all -operation g,f,f,f,f,f,f,f,f -nastring NA");
system("rm $ARGV[0].annovar.refGene* $ARGV[0].annovar.hg19*_dropped $ARGV[0].annovar.hg19*_filtered  $ARGV[0].annovar.*log $ARGV[0].annovar");

my %bamFile;
open(In, $ARGV[1]) or die;
while(<In>){
  chomp;
  my @array = split(/\t/, $_);
  # name, file
  $bamFile{$array[0]} = $array[1];
}
close(In);

# read depth
system("mkdir $ARGV[0].readN");
foreach my $name(sort keys %bamFile){
  print "$bamFile{$name}\n";
  system("/program/install/Anaconda3-2018.12/bin/python3.7 /WGS/script/count_depth_v2.py $bamFile{$name} $ARGV[0].bed > $ARGV[0].readN/$name.snv.readN");

  open(Out, ">$ARGV[0].readN/$name.snv");
  print Out "chr\tpos\tref_reads\tvar_reads\tvaf\n";
  open(In, "$ARGV[0].readN/$name.snv.readN");
  while(<In>){
    if(/^#/){ next; }
    chomp;
    my @array = split(/\t/, $_);
    my $tmpRefReadN = $array[5] - $array[6];
    my $tmpMutFreq = 0;
    if( $array[5]>0 ){
      $tmpMutFreq = $array[6]/$array[5]*100;
    }
    print Out "$array[0]\t$array[1]\t$tmpRefReadN\t$array[6]\t$tmpMutFreq\n";
  }
  close(In);
  close(Out);
}

