# cat  ../*CNA | awk -F '\t' '{print "Sample\t"  $1 "\t" $2 "\t" $3 "\tSample\t\t" $9}' - > BIC-Seq.SCNA

label="BIC-Seq"
path="/WGS/5-CopyNumberAlteration/"$label"/summary20240703/"
mkdir -p $path/CIN
cd $path/CIN

file_in=$path/CIN/BIC-Seq.SCNA
rm $file_in
for file in ../*.BICseq.SCNA; do
  # 提取文件名（去除路径，仅保留文件名本身）
  filename=$(basename "$file")
  # 处理当前文件，用文件名替换Sample
  awk -F '\t' -v fn="$filename" '{print fn "\t" $1 "\t" $2 "\t" $3 "\t" fn "\t\t" $9}' "$file" >> $file_in
done

file_len=/meta/ensemble_GRCh37.85/Homo_sapiens.GRCh37.85.protein_coding.len
file_size=/meta/hg19/Chromosomes.len
file_anno=ratio_${label}_anno.csv

sed '1d' $file_in|awk -F '\t' '{print $2"\t"$3"\t"$4"\t0\t-"}' >${label}.avinput
perl /meta/annovar/annotate_variation.pl -geneanno --neargene 2000 -buildver hg19 -dbtype refGene -outfile $label  -exonsort ${label}.avinput /meta/annovar/humandb
file1=${label}.variant_function
awk -F '\t' 'FILENAME==ARGV[1]{id=$3":"$4"-"$5;genes[id]=$2;regions[id]=$1}FILENAME==ARGV[2]{id1=$2":"$3"-"$4;print $2"\t"$3"\t"$4"\t"$5"\t"$7"\t"regions[id1]"\t"genes[id1]}' $file1 $file_in >$file_anno

sed '1d' $file_anno|awk -F '\t' '!($6~/downstream/)&&!($6~/upstream/)&&!($6~/UTR/)&&!($6~/intergenic/){print}'|awk -F '\t' '{split($7,a,",");for(i in a){print a[i]"\t"$4"\t"$5"\t"$1}}' >CNV_${label}_long_temp1.txt
awk -F '\t' 'FILENAME==ARGV[1]{lens[$1]=$2}FILENAME==ARGV[2]{print $1"\t"$2"\t"$3"\t"$4"\t"lens[$1]}' $file_len CNV_${label}_long_temp1.txt >CNV_${label}_long_temp.txt
awk -F '\t' 'FILENAME==ARGV[1]{sizes[$1]=$2}FILENAME==ARGV[2]{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"sizes[$4]}' $file_size CNV_${label}_long_temp.txt |awk -F '\t' '$4!=""{print }' >CNV_${label}_long_withChr.txt
cat CNV_${label}_long_withChr.txt|awk -F '\t' '{if($3<0){$3=-($3)};cin=$3*$5/$6;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"cin}'  >CNV_${label}_long_withCIN.txt
awk -F '\t' '{sums[$2":"$4]+=$7}END{for(i in sums){print i"\t"sums[i]} }' CNV_${label}_long_withCIN.txt|awk -F ':' '{print $1"\t"$2}'|sort -k1,1 -k2,2  >CNV_${label}_chr_CIN.txt
awk -F '\t' '$2!="chrX"&&$2!="chrY"&&$2!="chrM"{print }' CNV_${label}_chr_CIN.txt |awk -F '\t' '{cin[$1]+=$3}END{for(i in cin){print i"\t"cin[i]/22}}' >CNV_${label}_CIN.txt
