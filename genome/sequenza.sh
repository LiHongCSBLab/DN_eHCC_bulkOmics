case=$1
control=$2
HOME=/WGS/
REF=/meta/hg19/WholeGenomeFasta/genome.fa
GC=/meta/sequenza/hg19.gc50Base.wig.gz

cd $HOME/5-CopyNumberAlteration/sequenza
# Process a FASTA file to produce a GC Wiggle track file:
# /picb/bigdata/program/install/sequenza-utils/bin/sequenza-utils gc_wiggle −w 50 --fasta $REF -o $GC
# Process BAM and Wiggle files to produce a seqz file:
/picb/bigdata/program/install/sequenza-utils/bin/sequenza-utils bam2seqz -n $HOME/3-Preprocess/GATK/${control}.mkdup.realn.recal.bam -t $HOME/3-Preprocess/GATK/${case}.mkdup.realn.recal.bam --fasta $REF -gc $GC -o ${case}.seqz.gz
# Post-process by binning the original seqz file:
/picb/bigdata/program/install/sequenza-utils/bin/sequenza-utils seqz_binning --seqz ${case}.seqz.gz -w 50 -o ${case}.small.seqz.gz
# remove chrM
gunzip -c ${case}.small.seqz.gz| grep -v "chrM" > ${case}.small.v2.seqz
gzip ${case}.small.v2.seqz

/picb/bigdata/program/install/R-3.4.2/bin/Rscript $HOME/script/sequenza.R ${case}.small.v2.seqz.gz


