case=$1
control=$2
WGS_HOME=/WGS/

mkdir $WGS_HOME/5-CopyNumberAlteration/BIC-Seq/$case
cd $WGS_HOME/5-CopyNumberAlteration/BIC-Seq/$case

if [ ! -f $WGS_HOME/3-Preprocess/GATK/${case}.mkdup.realn.recal.bam  ]; then
  exit
fi
if [ ! -f $WGS_HOME/3-Preprocess/GATK/${control}.mkdup.realn.recal.bam  ]; then
  exit
fi

# get the mapping positions of the unique mapped reads
for ((i=1;i<=22;i++)); do
if [ ! -f ${case}.chr${i}.seq ]; then
  samtools view -q 10 -F 3332 -f 2 $WGS_HOME/3-Preprocess/GATK/${case}.mkdup.realn.recal.bam chr${i} | cut -f4 > ${case}.chr${i}.seq
fi
if [ ! -f ${control}.chr${i}.seq ]; then
  samtools view -q 10 -F 3332 -f 2 $WGS_HOME/3-Preprocess/GATK/${control}.mkdup.realn.recal.bam chr${i} | cut -f4 > ${control}.chr${i}.seq
fi  
done

perl $WGS_HOME/script/generate_BICseq_configure.pl $case $control .
if [ ! -f ${case}.chr${i}.norm.bin ]; then
  /picb/bigdata/program/install/NBICseq-norm_v0.2.4/NBICseq-norm.pl ${case}.BICseq-norm_configure ${case}.BICseq-norm.parameter
fi
if [ ! -f ${control}.chr${i}.norm.bin ]; then
  /picb/bigdata/program/install/NBICseq-norm_v0.2.4/NBICseq-norm.pl ${control}.BICseq-norm_configure ${control}.BICseq-norm.parameter
fi

# The log2 copy ratios for the tumor (normal) genome are defined as the log2 ratios between the observed and expected number of tumor (normal) reads in the segments. 
/picb/bigdata/program/install/NBICseq-seg_v0.7.2/NBICseq-seg.pl --fig=${case}.BICseq.SCNA.pdf --control --bootstrap ${case}.BICseq-seg_configure ${case}.BICseq.SCNA
# The CNV predictions are taken as regions with log2 tumor/expected ratio >0.2 or <−0.2.  Also filtered regions harboring likely germline variations, i.e. regions with log2 control/expected ratios >0.2 and <−0.2.

