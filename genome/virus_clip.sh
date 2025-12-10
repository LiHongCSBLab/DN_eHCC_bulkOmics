#!/bin/bash

# Usage: bash virus_clip_huisuan.sh $1 $2 $3

#tools
bwa=/program/bin/bwa
samtools=/program/bin/samtools
virus_clip=/program/install/Virus-Clip/virus_clip.pl
blastn=/program/install/ncbi-blast-2.2.29+/bin/blastn
annovar=/meta/annovar//annotate_variation.pl

#directory
SEQprefix=$1
OUTdir=$2

#resources
refFasta=/WGS/2-virusClip/HBV/HBV.fa
blastdb=/meta/hg19/WholeGenomeFasta/genome.fa
annovardb=/meta/annovar/humandb

#procedures
file=$3

mkdir ${OUTdir}/align
$bwa mem -t 2 ${refFasta} ${SEQprefix}.filter.R1.fastq ${SEQprefix}.filter.R2.fastq > ${OUTdir}/align/${file}.bwa.sam

#extract softclip reads
$samtools view -Sh ${OUTdir}/align/${file}.bwa.sam | awk '$6 ~ /S/ {print $2"\t"$4"\t"$5"\t"$6"\t"$10}' > ${OUTdir}/align/${file}.softclip.txt

#detect human and virus integration breakpoint and sequence
mkdir ${OUTdir}/result/
mkdir ${OUTdir}/result/${file}
cd ${OUTdir}/result/${file}
perl ${virus_clip} ${OUTdir}/align/${file}.softclip.txt ${blastn} ${blastdb} ${annovar} ${annovardb}

rm ${OUTdir}/align/${file}.bwa.sam

