#!/bin/sh
set -e
set -u

if [ $# -ne 1 ]; then
	echo "Usage: wgs.pipleline.hg19.sh NORMAL_ID"
	exit 1
fi

NORMAL_ID=$1

HOME=/HOME/
WGS_HOME=/WGS/
REF=${HOME}/meta/hg19/WholeGenomeFasta/genome.fa
bwaIndex=${HOME}/meta/hg19/BWAIndex/version0.7.12/genome.fa
LOG_HOME=${WGS_HOME}/log
if [ ! -d ${WGS_HOME}/log ]; then
        mkdir ${WGS_HOME}/log
fi
LOG=${LOG_HOME}/${NORMAL_ID}.log
tmp_HOME=${WGS_HOME}/tmp
if [ ! -d ${WGS_HOME}/tmp ]; then
        mkdir ${WGS_HOME}/tmp
fi

# QC
if [ ! -d ${WGS_HOME}/1-QC ]; then
	mkdir ${WGS_HOME}/1-QC 
fi
if [ ! -d ${WGS_HOME}/1-QC/FASTQC ]; then
        mkdir ${WGS_HOME}/1-QC/FASTQC
fi
cd ${WGS_HOME}/1-QC

echo "`date` : ${NORMAL_ID} start FASTQC " > ${LOG}
fastqc -o FASTQC ${WGS_HOME}/0-Raw/Sample_${NORMAL_ID}/*fastq.gz
echo "`date` : ${NORMAL_ID} finish FASTQC " >> ${LOG}

if [ ! -d filter ]; then
        mkdir filter
fi
if [ ! -d filter/${NORMAL_ID} ]; then
  mkdir filter/${NORMAL_ID}
fi
/home/lihong/script/fqclean -t 8 -o filter/$NORMAL_ID/$NORMAL_ID ${WGS_HOME}/0-Raw/Sample_${NORMAL_ID}/${NORMAL_ID}_combined_R1.fastq.gz ${WGS_HOME}/0-Raw/Sample_${NORMAL_ID}/${NORMAL_ID}_combined_R2.fastq.gz
cat filter/$NORMAL_ID/$NORMAL_ID*ovl.r1.fastq filter/$NORMAL_ID/$NORMAL_ID*gap.r1.fastq filter/$NORMAL_ID/$NORMAL_ID*.readthrough.r1.fastq > filter/$NORMAL_ID/${NORMAL_ID}.filter.R1.fastq
cat filter/$NORMAL_ID/$NORMAL_ID*ovl.r2.fastq filter/$NORMAL_ID/$NORMAL_ID*gap.r2.fastq filter/$NORMAL_ID/$NORMAL_ID*.readthrough.r2.fastq > filter/$NORMAL_ID/${NORMAL_ID}.filter.R2.fastq

# BWA
if [ ! -d ${WGS_HOME}/2-Mapping ]; then
	mkdir ${WGS_HOME}/2-Mapping
fi
cd ${WGS_HOME}/2-Mapping

if [ ! -d BWA ]; then
        mkdir BWA
fi

echo "`date` : ${NORMAL_ID} start bwa mem mapping" >> ${LOG}
bwa mem -t 16 $bwaIndex -r "@RG ID:${NORMAL_ID} LB:${NORMAL_ID} SM:${NORMAL_ID} PL:ILLUMINA PU:${NORMAL_ID}" ${WGS_HOME}/1-QC/filter/${NORMAL_ID}/${NORMAL_ID}.filter.R1.fastq ${WGS_HOME}/1-QC/filter/${NORMAL_ID}/${NORMAL_ID}.filter.R2.fastq > BWA/${NORMAL_ID}.sam
samtools view -bS BWA/${NORMAL_ID}.sam -o BWA/${NORMAL_ID}.bam 
samtools sort -m 3000M BWA/${NORMAL_ID}.bam -o BWA/${NORMAL_ID}.sort.bam
rm BWA/${NORMAL_ID}.sam BWA/${NORMAL_ID}.bam 
echo "`date` : ${NORMAL_ID} finish bwa mem mapping" >> ${LOG}

# Picard preprocessing for normal sample
if [ ! -d ${WGS_HOME}/3-Preprocess ]; then
        mkdir ${WGS_HOME}/3-Preprocess
fi
cd ${WGS_HOME}/3-Preprocess

if [ ! -d picard ]; then
        mkdir picard
fi

echo "`date` : ${NORMAL_ID} start marking duplicates" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME} -jar /picb/bigdata/program/install/picard-tools-2.3.0/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${WGS_HOME}/2-Mapping/BWA/${NORMAL_ID}.sort.bam OUTPUT=picard/${NORMAL_ID}.tmp.bam METRICS_FILE=picard/${NORMAL_ID}.mkdup.metrics MAX_RECORDS_IN_RAM=5000000
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME} -jar /picb/bigdata/program/install/picard-tools-2.3.0/picard.jar AddOrReplaceReadGroups I=picard/${NORMAL_ID}.tmp.bam O=picard/${NORMAL_ID}.mkdup.bam LB=${NORMAL_ID} PL=ILLUMINA PU=${NORMAL_ID} SM=${NORMAL_ID}
rm picard/${NORMAL_ID}.tmp.bam
echo "`date` : ${NORMAL_ID} finish marking duplicates" >> ${LOG}

echo "`date` : ${NORMAL_ID} start building sam index" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME} -jar /picb/bigdata/program/install/picard-tools-2.3.0/picard.jar BuildBamIndex VALIDATION_STRINGENCY=SILENT INPUT=picard/${NORMAL_ID}.mkdup.bam MAX_RECORDS_IN_RAM=5000000
echo "`date` : ${NORMAL_ID} finish building sam index" >> ${LOG}

echo "`date` : ${NORMAL_ID} start calculating mapping statistaics" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME} -jar /picb/bigdata/program/install/picard-tools-2.3.0/picard.jar CollectAlignmentSummaryMetrics R=$REF I=picard/${NORMAL_ID}.mkdup.bam O=picard/${NORMAL_ID}.mkdup.bam.summaryMetrics
echo "`date` : ${NORMAL_ID} finishing calculating mapping statistaics" >> ${LOG}

# very slow
# java -Xmx10g -Djava.io.tmpdir=${tmp_HOME} -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T DepthOfCoverage -R $REF -I picard/${NORMAL_ID}.mkdup.bam -o picard/${NORMAL_ID}.mkdup.bam.coverage

# samtools depth -a picard/${NORMAL_ID}.mkdup.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > picard/${NORMAL_ID}.mkdup.bam.depth

if [ ! -d GATK ]; then
        mkdir GATK
fi

echo "`date` : ${NORMAL_ID} start creating realignment target" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME}  -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -nt 8 -T RealignerTargetCreator -R $REF -I picard/${NORMAL_ID}.mkdup.bam  -o GATK/${NORMAL_ID}.mkdup.intervals --known ${HOME}/meta/GATK/1000G_phase1.indels.hg19.vcf --known ${HOME}/meta/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf
echo "`date` : ${NORMAL_ID} finish creating realignment target" >> ${LOG}

echo "`date` : ${NORMAL_ID} start realignment" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME}  -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I picard/${NORMAL_ID}.mkdup.bam -targetIntervals GATK/${NORMAL_ID}.mkdup.intervals -o GATK/${NORMAL_ID}.mkdup.realn.bam -known ${HOME}/meta/GATK/1000G_phase1.indels.hg19.vcf -known ${HOME}/meta/GATK/Mills_and_1000G_gold_standard.indels.hg19.vcf
echo "`date` : ${NORMAL_ID} finish realignment" >> ${LOG}

rm picard/${NORMAL_ID}.mkdup.bam picard/${NORMAL_ID}.mkdup.bai GATK/${NORMAL_ID}.mkdup.intervals


perl $WGS_HOME/script/removeG.pl GATK/${NORMAL_ID}

echo "`date` : ${NORMAL_ID} start creating base recalibrator" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME}  -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -nct 16 -T BaseRecalibrator -R $REF -I GATK/${NORMAL_ID}.mkdup.realn.removeG.bam -o GATK/${NORMAL_ID}.mkdup.realn.recal -knownSites ${HOME}/meta/GATK/dbsnp_138.hg19.vcf
echo "`date` : ${NORMAL_ID} finish creating base recalibrator" >> ${LOG}

echo "`date` : ${NORMAL_ID} start recalibration" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME}  -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -nct 16 -T PrintReads -R $REF -I GATK/${NORMAL_ID}.mkdup.realn.removeG.bam -o GATK/${NORMAL_ID}.mkdup.realn.recal.bam --BQSR GATK/${NORMAL_ID}.mkdup.realn.recal
echo "`date` : ${NORMAL_ID} finish recalibration" >> ${LOG}

echo "`date` : ${NORMAL_ID} start insertSize" >> ${LOG}
samtools view GATK/${NORMAL_ID}.mkdup.realn.recal.bam | awk '{print $9}' | grep -v "-" | sort -n | uniq -c > GATK/${NORMAL_ID}.mkdup.realn.recal.insertSize
echo "`date` : ${NORMAL_ID} finish insertSize" >> ${LOG}

rm GATK/${NORMAL_ID}.mkdup.realn.bam GATK/${NORMAL_ID}.mkdup.realn.bai GATK/${NORMAL_ID}.mkdup.realn.recal

# Germline SNP and INDEL calling
if [ ! -d ${WGS_HOME}/4-Germline ]; then
        mkdir ${WGS_HOME}/4-Germline
fi
cd ${WGS_HOME}/4-Germline

if [ ! -d GATK ]; then
        mkdir GATK
fi

echo "`date` : ${NORMAL_ID} start GATK calling" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME}  -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T UnifiedGenotyper -nct 2 -R $REF -glm BOTH --dbsnp ${HOME}/meta/GATK/dbsnp_138.hg19.vcf -I ${WGS_HOME}/3-Preprocess/GATK/${NORMAL_ID}.mkdup.realn.recal.bam -o GATK/${NORMAL_ID}_germline.GATK.raw.vcf
echo "`date` : ${NORMAL_ID} finish GATK calling" >> ${LOG}

echo "`date` : ${NORMAL_ID} start GATK SNP filtering" >> ${LOG}
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME} -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T SelectVariants -R $REF -selectType SNP --variant GATK/${NORMAL_ID}_germline.GATK.raw.vcf -o GATK/${NORMAL_ID}_germline.snp.GATK.raw.vcf
java -Xmx10g -Djava.io.tmpdir=${tmp_HOME} -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T VariantFiltration -R $REF --variant GATK/${NORMAL_ID}_germline.snp.GATK.raw.vcf -o GATK/${NORMAL_ID}_germline.snp.GATK.fil.vcf --filterName 'REJECT' --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
vcf-annotate --hard-filter GATK/${NORMAL_ID}_germline.snp.GATK.fil.vcf > GATK/${NORMAL_ID}_germline.snp.GATK.fil.passed.vcf
vcf-stats GATK/${NORMAL_ID}_germline.snp.GATK.fil.passed.vcf > GATK/${NORMAL_ID}_germline.snp.GATK.fil.passed.vcf.stat
java -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T VariantEval -R $REF -D ${HOME}/meta/GATK/dbsnp_138.hg19.vcf --eval GATK/${NORMAL_ID}_germline.snp.GATK.fil.passed.vcf -o GATK/${NORMAL_ID}_germline.snp.GATK.fil.passed.vcf.eval
echo "`date` : ${NORMAL_ID} finish GATK SNP filtering" >> ${LOG}

echo "`date` : ${NORMAL_ID} start GATK INDEL filtering" >> ${LOG}
java -Xmx10g -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T SelectVariants -R $REF -selectType INDEL --variant GATK/${NORMAL_ID}_germline.GATK.raw.vcf -o GATK/${NORMAL_ID}_germline.indel.GATK.raw.vcf
java -Xmx10g -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T VariantFiltration -R $REF --variant GATK/${NORMAL_ID}_germline.indel.GATK.raw.vcf -o GATK/${NORMAL_ID}_germline.indel.GATK.fil.vcf --filterName 'REJECT' --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
vcf-annotate --hard-filter GATK/${NORMAL_ID}_germline.indel.GATK.fil.vcf > GATK/${NORMAL_ID}_germline.indel.GATK.fil.passed.vcf
vcf-stats GATK/${NORMAL_ID}_germline.indel.GATK.fil.passed.vcf > GATK/${NORMAL_ID}_germline.indel.GATK.fil.passed.vcf.stat
java -jar /picb/bigdata/program/install/GATK-3.5/GenomeAnalysisTK.jar -T VariantEval -R $REF -D ${HOME}/meta/GATK/dbsnp_138.hg19.vcf --eval GATK/${NORMAL_ID}_germline.indel.GATK.fil.passed.vcf -o GATK/${NORMAL_ID}_germline.indel.GATK.fil.passed.vcf.eval
echo "`date` : ${NORMAL_ID} finish GATK INDEL filtering" >> ${LOG}

if [ ! -d ANNOVA ]; then
        mkdir ANNOVA
fi

echo "`date` : ${NORMAL_ID} start ANNOVAR SNP annotation" >> ${LOG}
perl ${HOME}/meta/annovar/convert2annovar.pl -format vcf4 GATK/${NORMAL_ID}_germline.snp.GATK.fil.passed.vcf -outfile ANNOVA/${NORMAL_ID}_germline.snp.GATK.fil.passed.annovar
perl ${HOME}/meta/annovar/table_annovar.pl  ANNOVA/${NORMAL_ID}_germline.snp.GATK.fil.passed.annovar ${HOME}/meta/annovar/humandb -buildver hg19  -protocol refGene,clinvar_20160302,snp138,cosmic70,esp6500siv2_all,ALL.sites.2014_10,exac03nontcga,ljb26_all -operation g,f,f,f,f,f,f,f -nastring NA
rm ANNOVA/${NORMAL_ID}_germline.snp.GATK.fil.annovar.refGene* ANNOVA/${NORMAL_ID}_germline.snp.GATK.fil.annovar.hg19*_dropped ANNOVA/${NORMAL_ID}_germline.snp.GATK.fil.annovar.hg19*_filtered
echo "`date` : ${NORMAL_ID} finish ANNOVAR SNP annotation" >> ${LOG}

echo "`date` : ${NORMAL_ID} start ANNOVAR INDEL annotation" >> ${LOG}
perl ${HOME}/meta/annovar/convert2annovar.pl -format vcf4 GATK/${NORMAL_ID}_germline.indel.GATK.fil.passed.vcf -outfile ANNOVA/${NORMAL_ID}_germline.indel.GATK.fil.passed.annovar
perl ${HOME}/meta/annovar/table_annovar.pl ANNOVA/${NORMAL_ID}_germline.indel.GATK.fil.passed.annovar ${HOME}/meta/annovar/humandb -buildver hg19  -protocol refGene,clinvar_20160302,snp138,cosmic70,esp6500siv2_all,ALL.sites.2014_10,exac03nontcga,ljb26_all -operation g,f,f,f,f,f,f,f  -nastring NA
rm ANNOVA/${NORMAL_ID}_germline.indel.GATK.fil.annovar.refGene* ANNOVA/${NORMAL_ID}_germline.indel.GATK.fil.annovar.hg19*_dropped ANNOVA/${NORMAL_ID}_germline.indel.GATK.fil.annovar.hg19*_filtered
echo "`date` : ${NORMAL_ID} finish ANNOVAR INDEL annotation" >> ${LOG}
