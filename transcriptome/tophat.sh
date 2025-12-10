#!/bin/bash

Seq_ID=$1
Sample_ID=$2
HOME=/rna_hg19/
GTF=/meta/ensemble_GRCh37.85/RSEM/Homo_sapiens.GRCh37.85.refined.gtf
bowtie2_index=//meta/ensemble_GRCh37.85/bowtie2/Homo_sapiens.GRCh37.dna


cd $HOME/1-QC/

if [ ! -d $HOME/1-QC/TrimGalore/ ]; then
  mkdir $HOME/1-QC/TrimGalore/
fi
/program/src/trim_galore_v0.4.4/trim_galore --phred33 --fastqc --illumina --output_dir TrimGalore --paired $HOME/0-Raw/Sample_${Seq_ID}/${Seq_ID}_combined_R1.fastq.gz $HOME/0-Raw/Sample_${Seq_ID}/${Seq_ID}_combined_R2.fastq.gz

cd $HOME/2-tophat-hg19
# pair-end
/program/bin/tophat2 -G $GTF -p 8 -o ${Sample_ID} ${bowtie2_index} $HOME/1-QC/TrimGalore/${Seq_ID}_combined_R1_val_1.fq.gz $HOME/1-QC/TrimGalore/${Seq_ID}_combined_R2_val_2.fq.gz

/program/bin/samtools index ${Sample_ID}/accepted_hits.bam

if [ ! -d $HOME/2-cufflinks-hg19 ]; then
  mkdir $HOME/2-cufflinks-hg19
fi
cd $HOME/2-cufflinks-hg19
if [ ! -e ${Sample_ID}/genes.fpkm_tracking ]; then
  /program/bin/cufflinks -G $GTF -p 8 -o ${Sample_ID} $HOME/2-tophat-hg19/${Sample_ID}/accepted_hits.bam
fi

#if [ ! -d $HOME/2-HTSeq-hg19 ]; then
#  mkdir $HOME/2-HTSeq-hg19
#fi
#cd $HOME/2-HTSeq-hg19
#samtools sort -n $HOME/2-tophat-hg19/${Sample_ID}/accepted_hits.bam -o $HOME/2-tophat-hg19/${Sample_ID}/accepted_hits_sorted.bam
#/program/anaconda3/bin/htseq-count -c ${Sample_ID}.count.csv  $HOME/2-tophat-hg19/${Sample_ID}/accepted_hits_sorted.bam $GTF

if [ ! -d $HOME/2-featureCounts-hg19 ]; then
  mkdir $HOME/2-featureCounts-hg19
fi
cd $HOME/2-featureCounts-hg19
/program/subread-2.0.6-Linux-x86_64/bin/featureCounts  -p --countReadPairs -t exon -g gene_id  -a  $GTF -o  ${Sample_ID}.count.csv  $HOME/2-tophat-hg19/${Sample_ID}/accepted_hits_sorted.bam



