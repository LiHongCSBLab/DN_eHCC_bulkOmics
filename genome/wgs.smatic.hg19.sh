#!/bin/sh
set -e
set -u

if [ $# -ne 2 ]; then
        echo "Usage: wgs.pipleline.hg19.sh NORMAL_ID TUMOR_ID"
        exit 1
fi

NORMAL_ID=$1
TUMOR_ID=$2
HOME=/HOME/
WGS_HOME=/WGS/
REF=${HOME}/meta/hg19/WholeGenomeFasta/genome.fa
LOG_HOME=${WGS_HOME}/log
LOG=${LOG_HOME}/${TUMOR_ID}_${NORMAL_ID}.somatic.pipeline.v0.log

if [ ! -d ${WGS_HOME}/5-SomaticMutation ]; then
        mkdir ${WGS_HOME}/5-SomaticMutation
fi

# Somatic SNV and INDEL calling by VarScan
cd ${WGS_HOME}/5-SomaticMutation/
if [ ! -d ${WGS_HOME}/5-SomaticMutation/varScan ]; then
        mkdir ${WGS_HOME}/5-SomaticMutation/varScan
fi
if [ ! -d ${WGS_HOME}/5-SomaticMutation/varScan_ANNOVAR ]; then
        mkdir ${WGS_HOME}/5-SomaticMutation/varScan_ANNOVAR
fi

echo "`date` : ${TUMOR_ID} start varScan calling" > ${LOG}
samtools mpileup -C50 -f $REF ${WGS_HOME}/3-Preprocess/GATK/${NORMAL_ID}.mkdup.realn.recal.bam ${WGS_HOME}/3-Preprocess/GATK/${TUMOR_ID}.mkdup.realn.recal.bam | awk '{if($4 >= 6) print $0}' | awk '{if($7 != 0) print $0}' | java -jar ${HOME}/program/bin/VarScan.v2.3.9.jar somatic --min-var-freq 0.05 --mpileup 1 --output-vcf 1  --output-snp varScan/${TUMOR_ID}_snv.varScan.vcf --output-indel varScan/${TUMOR_ID}_indel.varScan.vcf
echo "`date` : ${TUMOR_ID} finish varScan calling" >> ${LOG}

echo "`date` : ${TUMOR_ID} start varScan filtering" >> ${LOG}
java -jar ${HOME}/program/bin/VarScan.v2.3.9.jar somaticFilter varScan/${TUMOR_ID}_snv.varScan.vcf --min-strands2 2 --min-var-freq 0.05 --p-value 0.01 --indel-file varScan/${TUMOR_ID}_indel.varScan.vcf --output-file varScan/${TUMOR_ID}_snv.varScan.fil.passed.vcf
java -jar ${HOME}/program/bin/VarScan.v2.3.9.jar somaticFilter varScan/${TUMOR_ID}_indel.varScan.vcf --min-strands2 2 --min-var-freq 0.05 --p-value 0.01 --output-file varScan/${TUMOR_ID}_indel.varScan.fil.passed.vcf
echo "`date` : ${TUMOR_ID} finish varScan filtering" >> ${LOG}

echo "`date` : ${TUMOR_ID} start varScan extracting" >> ${LOG}
java -jar ${HOME}/program/bin/VarScan.v2.3.9.jar processSomatic varScan/${TUMOR_ID}_snv.varScan.fil.passed.vcf --min-tumor-freq 0.05 --p-value 0.01
java -jar ${HOME}/program/bin/VarScan.v2.3.9.jar processSomatic varScan/${TUMOR_ID}_indel.varScan.fil.passed.vcf --min-tumor-freq 0.05 --p-value 0.01
echo "`date` : ${TUMOR_ID} finish varScan extracting" >> ${LOG}

perl ${HOME}/meta/annovar/convert2annovar.pl -format vcf4old varScan/${TUMOR_ID}_snv.varScan.fil.passed.Somatic.hc.vcf -outfile varScan_ANNOVAR/${TUMOR_ID}_snv.varScan.fil.passed.Somatic.hc.annovar
perl ${HOME}/meta/annovar/table_annovar.pl  varScan_ANNOVAR/${TUMOR_ID}_snv.varScan.fil.passed.Somatic.hc.annovar ${HOME}/meta/annovar/humandb -buildver hg19  -protocol refGene,clinvar_20160302,snp138,cosmic70,esp6500siv2_all,ALL.sites.2014_10,exac03nontcga,ljb26_all -operation g,f,f,f,f,f,f,f -nastring NA
perl ${HOME}/meta/annovar/convert2annovar.pl -format vcf4old varScan/${TUMOR_ID}_indel.varScan.fil.passed.Somatic.hc.vcf -outfile varScan_ANNOVAR/${TUMOR_ID}_indel.varScan.fil.passed.Somatic.hc.annovar
perl ${HOME}/meta/annovar/table_annovar.pl  varScan_ANNOVAR/${TUMOR_ID}_indel.varScan.fil.passed.Somatic.hc.annovar ${HOME}/meta/annovar/humandb -buildver hg19  -protocol refGene,clinvar_20160302,snp138,cosmic70,esp6500siv2_all,ALL.sites.2014_10,exac03nontcga,ljb26_all -operation g,f,f,f,f,f,f,f -nastring NA

# Somatic SNV calling by muTect
cd ${WGS_HOME}/5-SomaticMutation/
if [ ! -d ${WGS_HOME}/5-SomaticMutation/muTect ]; then
        mkdir ${WGS_HOME}/5-SomaticMutation/muTect
fi
if [ ! -d ${WGS_HOME}/5-SomaticMutation/muTect_ANNOVAR ]; then
        mkdir ${WGS_HOME}/5-SomaticMutation/muTect_ANNOVAR
fi

echo "`date` : ${TUMOR_ID} start pre-processing BAM" >> ${LOG}
# ignore reads with too many mismatches or very low quality scores

echo "`date` : ${TUMOR_ID} finish pre-processing BAM" >> ${LOG}

echo "`date` : ${TUMOR_ID} start muTect calling" >> ${LOG}
/usr/bin/java -Xmx10g -jar ${HOME}/program/install/mutect-1.1.7.jar --analysis_type MuTect --reference_sequence $REF --cosmic ${HOME}/meta/GATK/cosmic_v74_20151029.hg19.vcf --dbsnp ${HOME}/meta/GATK/dbsnp_138.hg19.vcf --input_file:normal ${WGS_HOME}/3-Preprocess/GATK/${NORMAL_ID}.mkdup.realn.recal.bam --input_file:tumor ${WGS_HOME}/3-Preprocess/GATK/${TUMOR_ID}.mkdup.realn.recal.bam --normal_sample_name NORMAL --tumor_sample_name TUMOR --out muTect/${TUMOR_ID}_snv.muTect.out --vcf muTect/${TUMOR_ID}_snv.muTect.vcf --fraction_contamination 0.2 --min_qscore 20 --gap_events_threshold 3
echo "`date` : ${TUMOR_ID} finish muTect calling" >> ${LOG}

echo "`date` : ${TUMOR_ID} start muTect filtering" >> ${LOG}
grep -v "REJECT" muTect/${TUMOR_ID}_snv.muTect.vcf > muTect/${TUMOR_ID}_snv.muTect.hc.vcf
echo "`date` : ${TUMOR_ID} finish muTect filtering" >> ${LOG}

perl ${HOME}/meta/annovar/convert2annovar.pl -format vcf4old muTect/${TUMOR_ID}_snv.muTect.hc.vcf -outfile muTect_ANNOVAR/${TUMOR_ID}_snv.muTect.hc.annovar
perl ${HOME}/meta/annovar/table_annovar.pl  muTect_ANNOVAR/${TUMOR_ID}_snv.muTect.hc.annovar ${HOME}/meta/annovar/humandb -buildver hg19  -protocol refGene,clinvar_20160302,snp138,cosmic70,esp6500siv2_all,ALL.sites.2014_10,exac03nontcga,ljb26_all -operation g,f,f,f,f,f,f,f -nastring NA

exit

## Somatic SNV and indel calling by Strelka
cd ${WGS_HOME}/5-SomaticMutation/
if [ ! -d ${WGS_HOME}/5-SomaticMutation/Strelka ]; then
        mkdir ${WGS_HOME}/5-SomaticMutation/Strelka
fi
#perl ${HOME}/program/install/strelka_workflow-1.0.14/bin/configureStrelkaWorkflow.pl --tumor=${WGS_HOME}/3-Preprocess/GATK/${TUMOR_ID}.mkdup.realn.recal.bam --normal=${WGS_HOME}/3-Preprocess/GATK/${NORMAL_ID}.mkdup.realn.recal.bam --ref=$REF --config=${HOME}/program/install/strelka_workflow-1.0.14/demo/strelka_demo_config.ini --output-dir=Strelka/${TUMOR_ID}


/picb/bigdata//program/install/strelka-2.7.1.centos5_x86_64/bin//configureStrelkaSomaticWorkflow.py --tumorBam=${WGS_HOME}/3-Preprocess/GATK/${TUMOR_ID}.mkdup.realn.recal.bam --normalBam=${WGS_HOME}/3-Preprocess/GATK/${NORMAL_ID}.mkdup.realn.recal.bam --referenceFasta=$REF --runDir=Strelka/${TUMOR_ID}
Strelka/${TUMOR_ID}/runWorkflow.py -m local -j 1 -g 4


## Somatic SNV calling by SomaticSniper
cd ${WGS_HOME}/5-SomaticMutation/
if [ ! -d ${WGS_HOME}/5-SomaticMutation/SomaticSniper ]; then
        mkdir ${WGS_HOME}/5-SomaticMutation/SomaticSniper
fi
# bam-somaticsniper [options] -f <ref.fasta> <tumor.bam> <normal.bam> <snv_output_file>
${HOME}/program/install/somatic-sniper-1.0.4/build/bin/bam-somaticsniper -q 1 -f $REF -F vcf ${WGS_HOME}/3-Preprocess/GATK/${TUMOR_ID}.mkdup.realn.recal.bam ${WGS_HOME}/3-Preprocess/GATK/${NORMAL_ID}.mkdup.realn.recal.bam SomaticSniper/${TUMOR_ID}.somaticsniper
# generate a samtools pileup (not mpileup) indel file
samtools mpileup -C50 -f $REF ${WGS_HOME}/3-Preprocess/GATK/${TUMOR_ID}.mkdup.realn.recal.bam > SomaticSniper/${TUMOR_ID}.indel_pileup
# Filter on standard filters using the indel file 
perl ${HOME}/program/install/somatic-sniper-1.0.4/src/scripts/snpfilter.pl -snp-file SomaticSniper/${TUMOR_ID}.somaticsniper -indel-file SomaticSniper/${TUMOR_ID}.indel_pileup
# Adapt the remainder for use with bam-readcount
perl ${HOME}/program/install/somatic-sniper-1.0.4/src/scripts/prepare_for_readcount.pl -snp-file SomaticSniper/${TUMOR_ID}.somaticsniper.SNPfilter
# Run bam-readcount
${HOME}/program/install/bam-readcount-master/build/bin/bam-readcount -q 1 -b 15 -f $REF -l SomaticSniper/${TUMOR_ID}.somaticsniper.SNPfilter.pos ${WGS_HOME}/3-Preprocess/GATK/${TUMOR_ID}.mkdup.realn.recal.bam > SomaticSniper/${TUMOR_ID}.readcounts.rc
# false positive filter
perl ${HOME}/program/install/somatic-sniper-1.0.4/src/scripts/fpfilter.pl -snp-file SomaticSniper/${TUMOR_ID}.somaticsniper.SNPfilter -readcount-file SomaticSniper/${TUMOR_ID}.readcounts.rc
# run the "high confidence" filter which filters based on the Somatic Score and mapping quality 
perl ${HOME}/program/install/somatic-sniper-1.0.4/src/scripts/highconfidence.pl -snp-file SomaticSniper/${TUMOR_ID}.somaticsniper.SNPfilter.fp_pass

## Somatic SNV calling by RADIA (RNA and DNA integrated analysis)

