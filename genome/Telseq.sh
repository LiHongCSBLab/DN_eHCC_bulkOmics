tumor=$1
WGS_HOME=/WGS

# The software is implemented in C++ (GCC version 4.8 or above)

/program/src/telseq-master/src/Telseq/telseq $WGS_HOME/3-Preprocess/GATK/${tumor}.mkdup.realn.recal.bam -o $WGS_HOME/5-telomere/Telseq/${tumor}.Telseq


