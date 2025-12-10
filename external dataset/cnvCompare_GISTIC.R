library(readxl)
library(ggplot2)
library(ggsci)

plot_GISTIC <- function( scores, title ){
  df = read.table("../hg19_Chromosomes.len", header=F, sep="\t")
  colnames(df) = c("chromName", "chromlength")
  df = df[-c(1,24:25), ]
  df$chromNum <- 1:length(df$chromName) 
  df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度
  # 得到每条染色体从0开始的起始坐标
  df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])
  # 计算每条染色体中间位置坐标，用来最后加文字
  tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
  df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle

  chromID <- scores$Chromosome
  scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
  scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
  scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1

  ggplot(scores, aes(StartPos, frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T))+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=1, color="gray")+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=0.6,label=sub("chr","",chromName)), size=3)+
  theme_minimal() + theme(axis.title.x = element_blank(),  axis.text.x = element_blank() ) + 
  ggtitle(title)
}

data1_seg = as.data.frame(read_excel("../38355797/Copy_Number_Alteration_20250815.xlsx"))
tmp = data.frame(data1_seg[,1:2], as.numeric(data1_seg[,3]), as.numeric(data1_seg[,4]), as.numeric(data1_seg[,4])-as.numeric(data1_seg[,3]), log2(as.numeric(data1_seg[,5]))-1 )
colnames(tmp) = c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean")
write.table( tmp, "Copy_Number_Alteration_20250815.seg", sep="\t", quote=F, row.names=F )

setwd("-ta 0.3 -td 0.3 -conf 0.9")

pdf("cnvCompare_GISTIC.pdf", width=10, height=3)
data0_scores = read.csv("NinN_HCC_GISTIC/scores.gistic", header=T, sep="\t")
plot_GISTIC( data0_scores, "NinN_HCC" )
data1_scores = read.csv("Copy_Number_Alteration_20250815_GISTIC/scores.gistic", header=T, sep="\t")
plot_GISTIC( data1_scores, "Copy_Number_Alteration_20250815" )
data2_scores = read.csv("TCGA_cna_hg19_GISTIC/scores.gistic", header=T, sep="\t")
plot_GISTIC( data2_scores, "TCGA_cna_hg19" )
dev.off()

library(VennDiagram)
amp_conf90 = list()
amp_conf90[["NinN_HCC"]] = unlist(strsplit(readLines("NinN_HCC_GISTIC/amp_genes.conf_90.txt")[1], split="\t"))[-1]
amp_conf90[["Copy_Number_Alteration_20250815"]] = unlist(strsplit(readLines("Copy_Number_Alteration_20250815_GISTIC/amp_genes.conf_90.txt")[1], split="\t"))[-1]
amp_conf90[["TCGA_cna_hg19"]] = unlist(strsplit(readLines("TCGA_cna_hg19_GISTIC/amp_genes.conf_90.txt")[1], split="\t"))[-1]
venn.diagram( amp_conf90, filename="cnvCompare_GISTIC_amp90.tiff" )
venn.diagram( sapply(amp_conf90, function(x){ unique(sub("\\..*", "", x)) }), filename="cnvCompare_GISTIC_amp90_v2.tiff" )

del_conf90 = list()
del_conf90[["NinN_HCC"]] = unlist(strsplit(readLines("NinN_HCC_GISTIC/del_genes.conf_90.txt")[1], split="\t"))[-1]
del_conf90[["Copy_Number_Alteration_20250815"]] = unlist(strsplit(readLines("Copy_Number_Alteration_20250815_GISTIC/del_genes.conf_90.txt")[1], split="\t"))[-1]
del_conf90[["TCGA_cna_hg19"]] = unlist(strsplit(readLines("TCGA_cna_hg19_GISTIC/del_genes.conf_90.txt")[1], split="\t"))[-1]
venn.diagram( del_conf90, filename="cnvCompare_GISTIC_del90.tiff" )
venn.diagram( sapply(del_conf90, function(x){ unique(sub("\\..*", "", x)) }), filename="cnvCompare_GISTIC_del90_v2.tiff" )

amp_broad = matrix(0, nrow=44, ncol=3)
colnames(amp_broad) = c("NinN_HCC", "CLCA", "TCGA")
rownames(amp_broad) = paste(rep(1:22,each=2), c("p", "q"), sep="")
del_broad = amp_broad
tmp = as.matrix(read.csv("NinN_HCC_GISTIC/broad_significance_results.txt", header=T, sep="\t", row.names=1))
amp_broad[rownames(tmp),"NinN_HCC"] = tmp[, "Amp.frequency"]
del_broad[rownames(tmp),"NinN_HCC"] = tmp[, "Del.frequency"]
tmp = as.matrix(read.csv("Copy_Number_Alteration_20250815_GISTIC/broad_significance_results.txt", header=T, sep="\t", row.names=1))
amp_broad[rownames(tmp),"CLCA"] = tmp[, "Amp.frequency"]
del_broad[rownames(tmp),"CLCA"] = tmp[, "Del.frequency"]
tmp = as.matrix(read.csv("TCGA_cna_hg19_GISTIC/broad_significance_results.txt", header=T, sep="\t", row.names=1))
amp_broad[rownames(tmp),"TCGA"] = tmp[, "Amp.frequency"]
del_broad[rownames(tmp),"TCGA"] = tmp[, "Del.frequency"]
library(pheatmap)
pdf("cnvCompare_broad.pdf")
pheatmap(amp_broad, cluster_cols=F, scale="none", main="amp_broad")
write.csv(amp_broad, "GISTIC_AMP_broad.csv")
pheatmap(del_broad, cluster_cols=F, scale="none", main="del_broad")
write.csv(del_broad, "GISTIC_DEL_broad.csv")

plot(1:4, rep(0,4), xaxt="n", ylim=c(0, 0.6), xlab="", ylab="Amp.frequency")
for( i in rownames(amp_broad) ){
  lines(1:3, amp_broad[i, ]  )
  text(3.5, amp_broad[i, 3], paste(i, amp_broad[i, 3], sep=" "), cex=0.8 )
}
axis(1, 1:3, labels=colnames(amp_broad) )

plot(1:4, rep(0,4), xaxt="n", ylim=c(0, 0.6), xlab="", ylab="Amp.frequency (>0.2 in CLCA and TCGA)")
for( i in rownames(amp_broad)[amp_broad[,2]>0.2 & amp_broad[,3]>0.2] ){
  lines(1:3, amp_broad[i, ]  )
  text(3.5, amp_broad[i, 3], paste(i, amp_broad[i, 3], sep=" "), cex=0.8 )
}
axis(1, 1:3, labels=colnames(amp_broad) )

plot(1:4, rep(0,4), xaxt="n", ylim=c(0, 0.6), xlab="", ylab="Del.frequency")
for( i in rownames(del_broad) ){
  lines(1:3, del_broad[i, ]  )
  text(3.5, del_broad[i, 3], paste(i, del_broad[i, 3], sep=" "), cex=0.8 )
}
axis(1, 1:3, labels=colnames(del_broad) )

plot(1:4, rep(0,4), xaxt="n", ylim=c(0, 0.6), xlab="", ylab="Del.frequency (>0.2 in CLCA and TCGA)")
for( i in rownames(del_broad)[del_broad[,2]>0.2 & del_broad[,3]>0.2] ){
  lines(1:3, del_broad[i, ]  )
  text(3.5, del_broad[i, 3], paste(i, del_broad[i, 3], sep=" "), cex=0.8 )
}
axis(1, 1:3, labels=colnames(del_broad) )

dev.off()

save.image("cnvCompare_GISTIC.Rdata")







