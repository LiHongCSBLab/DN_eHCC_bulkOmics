#  cut -f3,5,6,4,13 /picb/bigdata/meta/annovar/humandb/hg19_refGene.txt | sort|uniq > exonic.bed
#
my @files = (`ls /picb/bigdata/project/HuiLiJian_DN_2022/WGS/5-SomaticMutation/muTect/*muTect.hc.vcf.gz`, `ls /picb/bigdata/project/HuiLiJian_DN_2022/WGS/5-SomaticMutation/varScan/*_snv.varScan.fil.passed.Somatic.hc.vcf.gz`);
# @files = (@files, `ls /picb/bigdata/project/Hui_cellline/DN/WGS/5-SomaticMutation/muTect/WGC045*muTect.hc.vcf.gz`, `ls /picb/bigdata/project/Hui_cellline/DN/WGS/5-SomaticMutation/varScan/WGC045*Somatic.hc.vcf.gz`);
# @files = (`ls /picb/bigdata/project/Hui_cellline/DN/WGS/5-SomaticMutation/varScan/*17s22033*_snv.varScan.fil.passed.Somatic.vcf.gz`);

chomp(@files);
foreach my $file(@files){
  my $name = $file;
  $name=~s/.*\///;
  $name=~s/\.vcf\.gz//;
  system("vcftools --gzvcf $file --bed exonic.bed --recode --out exonic_vcf/$name");
}
system("cat exonic_vcf/*recode.vcf|grep -v '#' |cut -f1,2,4,5|sort|uniq > exonic_mut.bed");
my $cmd = "awk '{print \$1\"\\t\"\$2\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}' exonic_mut.bed > exonic_mut.annovar" ;
system($cmd);

system("perl /picb/bigdata/meta/annovar/table_annovar.pl exonic_mut.annovar /picb/bigdata/meta/annovar/humandb -buildver hg19 -protocol refGene,clinvar_20160302,snp138,cosmic70,esp6500siv2_all,ALL.sites.2014_10,exac03nontcga,ljb26_all -operation g,f,f,f,f,f,f,f -nastring NA");
system("rm exonic_mut.annovar*_dropped exonic_mut.annovar*_filtered exonic_mut.annovar*refGene*");
my $cmd = "grep -P \"\texonic\t\" exonic_mut.annovar.hg19_multianno.txt | grep -P -v \"\\tsynonymous\" > exonic_functionalMut.annovar.hg19_multianno.txt";
system($cmd);
my $cmd = "grep -P \"\tsplicing\t\" exonic_mut.annovar.hg19_multianno.txt >> exonic_functionalMut.annovar.hg19_multianno.txt";
system($cmd);
system("cut -f1,2,4,5 exonic_functionalMut.annovar.hg19_multianno.txt > exonic_functionalMut.bed ");

my @bamFiles = `ls /picb/bigdata/project/HuiLiJian_DN_2022/WGS/3-Preprocess/GATK/*realn.recal.bam`;
chomp(@bamFiles);
foreach my $bamFile (@bamFiles){
  if( $bamFile=~/R19059763LD01-18s18941-eHCC1_20200306/ ){ next; }
  if( $bamFile=~/DN1705090010/ ){ next; }
  my $name = $bamFile;
  $name=~s/.*\///;
  $name=~s/\.mkdup\.realn\.recal\.bam//;
  system("/picb/bigdata/program/install/Anaconda3-2018.12/bin/python3.7 /picb/bigdata/project/HuiLiJian_DN_2022/WGS/script/count_depth_v2.py $bamFile exonic_functionalMut.bed > exonic_mut_readN/$name.readN");
}

my %readN; my @names;
my @files = `ls exonic_mut_readN/*.readN`;
chomp(@files);
foreach my $file(@files){
  $name = $file;
  $name =~s/.*\///;
  $name=~s/\.readN//;
  push @names, $name;
  open(In, $file);
  while(<In>){
    chomp;
    if(/^#/){ next; }
    my @array = split(/\t/, $_);
    $readN{"$array[0]_$array[1]_$array[2]_$array[3]"}{$name} = "$array[6]/$array[5]";
  }
  close(In);
}

open(In, "exonic_functionalMut.annovar.hg19_multianno.txt");
open(Out, ">exonic_functionalMut.readN");
print Out "Chr\tStart\tRef\tAlt\tGene\tExonicFunc\tsnp138\tcosmic70\tesp6500siv2_all\tALL.sites.2014_10\tExAC_nontcga_ALL";
foreach my $name(@names){
  print Out "\t$name";
}
print Out "\n";
while(<In>){
  chomp;
  my @array = split(/\t/, $_);
  print Out  "$array[0]\t$array[1]\t$array[3]\t$array[4]\t$array[6]\t$array[8]\t$array[15]\t$array[16]\t$array[17]\t$array[18]\t$array[19]";
  foreach my $name(@names){
    print Out "\t".$readN{"$array[0]_$array[1]_$array[3]_$array[4]"}{$name};
  }
  print Out "\n";
}
close(In);
close(Out);

