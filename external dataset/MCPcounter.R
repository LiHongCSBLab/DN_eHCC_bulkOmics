################### TCGA FPKM to TPM #############
file = "TCGA-LIHC-FPKM.csv"
data = read.csv(file, header=T, row.names=1)

library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_transcript_info <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "transcript_length"),  
  filters = "ensembl_gene_id",
  values = rownames(data),
  mart = mart
)
gene_transcript_info = gene_transcript_info[gene_transcript_info[, "gene_biotype"]=="protein_coding", ]

# 按 Gene ID 取最长转录本长度（作为基因长度近似值）
tmp <- aggregate(
  transcript_length ~ ensembl_gene_id,
  data = gene_transcript_info,
  FUN = max  # 取最长转录本（最接近外显子总长度）
)
gene_length = tmp[,2]
names(gene_length) = tmp[,1]

fpkm_to_tpm <- function(fpkm_matrix, gene_length) {
  # 验证输入：基因ID匹配、无负数值
  stopifnot(all(rownames(fpkm_matrix) == names(gene_length)))  # 基因ID必须一致
  if (any(fpkm_matrix < 0)) fpkm_matrix[fpkm_matrix < 0] <- 0  # 过滤负数值（测序误差）
  
  # 步骤1：计算每个基因的 FPKM / 基因长度（kb）→ 即 FPKM/(L/1000)
  fpkm_per_kb <- fpkm_matrix / (gene_length / 1000)  # 每行对应一个基因的长度校正值
  
  # 步骤2：计算每个样本的缩放因子（总和/1e6）
  scaling_factor <- colSums(fpkm_per_kb) / 1e6
  
  # 步骤3：计算TPM（避免除以0）
  tpm_matrix <- fpkm_per_kb / scaling_factor
  return(tpm_matrix)
}

tpm <- fpkm_to_tpm(data[names(gene_length), ], gene_length)
tmp = unique(gene_transcript_info[,1:3])
tmp1 = tmp[,3]
names(tmp1) = tmp[,1]
tmp2 = tmp1[rownames(tpm)]
tpm_v2 = tpm[tmp2!="" & tmp2!="PINX1", ]
rownames(tpm_v2) = tmp2[tmp2!="" & tmp2!="PINX1"]

colnames(tpm_v2) = gsub("[01]1A.*", "01", colnames(tpm_v2))
colnames(tpm_v2) = gsub("[01]1B.*", "01", colnames(tpm_v2))
samples = colnames(read.table("cBioPortal_TCGA_mrna_seq_v2_rsem.txt", sep="\t", header=T))[-(1:2)]
tpm_v2 = tpm_v2[, intersect(samples, colnames(tpm_v2))]
write.csv(tpm_v2, "TCGA-LIHC-TPM.csv")

############## run MCPcounter ###############
library(MCPcounter)
genesM = read.table("MCPcounter_SignatureGenes.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)

file = "TCGA-LIHC-TPM.csv"
data = read.csv(file, header=T, row.names=1)
res = MCPcounter.estimate(data, featuresType=c("HUGO_symbols"), genes=genesM )
write.table(res, paste(file, ".MCPcounter", sep=""), quote=F, sep="\t", row.names=T, col.names=T)

file = "94Sample.ProteinCoding.tpm.csv"
tmp = read.csv(file, sep=",", header=T)
tmp = tmp[ tapply(1:nrow(tmp), tmp[,3], function(x){ if( length(x)>1 ){ x[1] }else{ x } }) , ]
data = tmp[,-(1:3)] 
rownames(data) = tmp[,3]
res = MCPcounter.estimate(data, featuresType=c("HUGO_symbols"), genes=genesM )
write.table(res, paste(file, ".MCPcounter", sep=""), quote=F, sep="\t", row.names=T, col.names=T)

file = "cBioPortal_CLCA_mrna_seq_tpm.txt"
data = read.table(file,sep="\t", header=T, row.names=1)
res = MCPcounter.estimate(data, featuresType=c("HUGO_symbols"), genes=genesM )
write.table(res, paste(file, ".MCPcounter", sep=""), quote=F, sep="\t", row.names=T, col.names=T)

file = "cBioPortal_CLCA_mrna_seq_tpm_zscores_ref_all_samples.txt"
data = read.table(file,sep="\t", header=T, row.names=1)
res = MCPcounter.estimate(data, featuresType=c("HUGO_symbols"), genes=genesM )
write.table(res, paste(file, ".MCPcounter", sep=""), quote=F, sep="\t", row.names=T, col.names=T)

file = "cBioPortal_TCGA_mrna_seq_v2_rsem.txt"
tmp = read.table(file,sep="\t", header=T)
tmp = tmp[tmp[,1]!="", ]
tmp = tmp[ tapply(1:nrow(tmp), tmp[,1], function(x){ if( length(x)>1 ){ x[1] }else{ x } }) , ]
data = as.matrix(tmp[,-(1:2)])
rownames(data) = tmp[,1]
res = MCPcounter.estimate(data, featuresType=c("HUGO_symbols"), genes=genesM )
write.table(res, paste(file, ".MCPcounter", sep=""), quote=F, sep="\t", row.names=T, col.names=T)

file = "cBioPortal_TCGA_mrna_seq_v2_rsem_zscores_ref_all_samples.txt"
tmp = read.table(file,sep="\t", header=T)
tmp = tmp[tmp[,1]!="", ]
tmp = tmp[ tapply(1:nrow(tmp), tmp[,1], function(x){ if( length(x)>1 ){ x[1] }else{ x } }) , ]
data = as.matrix(tmp[,-(1:2)])
rownames(data) = tmp[,1]
res = MCPcounter.estimate(data, featuresType=c("HUGO_symbols"), genes=genesM )
write.table(res, paste(file, ".MCPcounter", sep=""), quote=F, sep="\t", row.names=T, col.names=T)

##################### compare results ###############
data0 = as.matrix(read.csv("94Sample.ProteinCoding.tpm.csv.MCPcounter", sep="\t", header=T))
data1 = as.matrix(read.csv("cBioPortal_CLCA_mrna_seq_tpm.txt.MCPcounter", sep="\t", header=T))
data2 = as.matrix(read.csv("TCGA-LIHC-TPM.csv.MCPcounter", sep="\t", header=T))
pdf("MCPcounter.pdf")
for(i in rownames(data0)){
  boxplot( list(NT=data0[i, colnames(data0)[grep("_NT", colnames(data0))]], DN=data0[i, colnames(data0)[grep("_DN", colnames(data0))]], eHCC=data0[i, colnames(data0)[grep("_T", colnames(data0))]], CLCA=data1[i, ], TCGA=data2[i, ]), ylab=i )
}
dev.off()

####################################
tmp = read.csv("94Sample.ProteinCoding.tpm.csv", sep=",", header=T)
tmp = tmp[ tapply(1:nrow(tmp), tmp[,3], function(x){ if( length(x)>1 ){ x[1] }else{ x } }) , ]
data0 = tmp[,-(1:3)] 
rownames(data0) = tmp[,3]
data1 = as.matrix(read.csv("cBioPortal_CLCA_mrna_seq_tpm.txt", sep="\t", header=T, row.names=1))
data2 = as.matrix(read.csv("TCGA-LIHC-TPM-424.csv", sep=",", header=T, row.names=1))
library(sva)
genes = intersect(intersect(rownames(data0), rownames(data1)), rownames(data2))
data_combined = cbind(data0[genes,], data1[genes,], data2[genes,])
batch <- factor(c(rep("data0", ncol(data0)), rep("data1", ncol(data1)), rep("data2", ncol(data2)))) 
data_corrected <- ComBat(dat = data_combined, batch = batch)
res_combined = MCPcounter.estimate(data_corrected, featuresType=c("HUGO_symbols"), genes=genesM )
write.table(res_combined, "RNAseqCorrected.MCPcounter", quote=F, sep="\t", row.names=T, col.names=T)

tmp = as.matrix(read.csv("sampleInfor", header=T, sep="\t"))
ids = tmp[,"finalID"]
ids = gsub("_Ca", "_T", ids)
ids = c(ids, paste( unique(gsub("_.*", "", ids)), "_NT", sep="") )

pdf("Corrected_MCPcounter.pdf", width=9)
p = matrix(nrow=nrow(res_combined), ncol=6)
rownames(p) = rownames(res_combined)
colnames(p) = c("out.VS.CLCA", "our.VS.TCGA", "CLCA.VS.TCGA", "NT.VS.DN", "TCGA_NT.VS.DN", "DN.VS.eHCC")
for(i in rownames(res_combined)){
  boxplot( list(NT=res_combined[i, intersect(ids, colnames(data0)[grep("_NT", colnames(data0))])], 
    TCGA_NT=res_combined[i, colnames(data2)[grep("\\.11A", colnames(data2))]], 
    DN=res_combined[i, intersect(ids, colnames(data0)[grep("_DN", colnames(data0))])], 
    eHCC=res_combined[i, intersect(ids, colnames(data0)[grep("_T", colnames(data0))])], 
    CLCA=res_combined[i, colnames(data1)], 
    TCGA_HCC=res_combined[i, colnames(data2)[grep("\\.01A", colnames(data2))]] ), 
    ylab=i )

  p[i, ] = c(
    wilcox.test( res_combined[i, intersect(ids, colnames(data0)[grep("_T", colnames(data0))])], res_combined[i, colnames(data1)] )$p.value,
    wilcox.test( res_combined[i, intersect(ids, colnames(data0)[grep("_T", colnames(data0))])], res_combined[i, colnames(data2)[grep("\\.01A", colnames(data2))]]  )$p.value,
    wilcox.test( res_combined[i, colnames(data1)], res_combined[i, colnames(data2)[grep("\\.01A", colnames(data2))]]  )$p.value,
    wilcox.test( res_combined[i, intersect(ids, colnames(data0)[grep("_NT", colnames(data0))])], res_combined[i, intersect(ids, colnames(data0)[grep("_DN", colnames(data0))])] )$p.value,
    wilcox.test( res_combined[i, colnames(data2)[grep("\\.11A", colnames(data2))]], res_combined[i, intersect(ids, colnames(data0)[grep("_DN", colnames(data0))])] )$p.value,
    wilcox.test( res_combined[i, intersect(ids, colnames(data0)[grep("_DN", colnames(data0))])] , res_combined[i, intersect(ids, colnames(data0)[grep("_T", colnames(data0))])] )$p.value
  )
}
dev.off()
write.csv(p, "Corrected_MCPcounter_P.csv")


cells = c("T cells", "CD8 T cells", "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Neutrophils", "Endothelial cells", "Fibroblasts")
library(reshape2)
df = rbind(
data.frame( melt( res_combined[cells, intersect(ids, colnames(data0)[grep("_NT", colnames(data0))]) ] ), dataset="our_NT" ),
data.frame( melt( res_combined[cells, intersect(ids, colnames(data0)[grep("_DN", colnames(data0))]) ] ), dataset="our_DN" ),
data.frame( melt( res_combined[cells, intersect(ids, colnames(data0)[grep("_T", colnames(data0))]) ] ), dataset="our_HCC" ),
data.frame( melt( res_combined[cells, colnames(data1) ] ), dataset="CLCA_HCC" ),
data.frame( melt( res_combined[cells, colnames(data2)[grep("\\.01A", colnames(data2))] ] ), dataset="TCGA_HCC" ) )
colnames(df) = c("Cell", "Sample", "Percentage", "Dataset")
df$Dataset <- factor(df$Dataset, levels = c("our_NT", "our_DN", "our_HCC", "CLCA_HCC", "TCGA_HCC"))
tmp <- tapply(1:nrow(df), df[, "Cell"], function(x){ y=df[x, "Percentage"]; (y-mean(y))/sd(y) })
tmpI <- tapply(1:nrow(df), df[, "Cell"], function(x){ x })
df[unlist(tmpI), "Percentage"] <- unlist(tmp)
pdf("Corrected_MCPcounter_boxplot.pdf", width=12, height=5)
library(ggplot2)
ggplot(df, aes(x = Cell, y = Percentage, fill = Dataset)) +
  # 箱线图（优先绘制，散点叠在上方）
  geom_boxplot(
    position = position_dodge(0.8),
    width = 0.7,
    outlier.shape = NA,  
    outlier.alpha = 0.5,  # 异常值透明度
    notch = TRUE,         # 凹槽（展示中位数95%置信区间）
    notchwidth = 0.8      # 凹槽宽度
  ) + 
  # 基础样式
  labs(
    x = "Cell",
    y = "Normalized proportion of cells",
    fill = "Dataset" # 图例标题
  ) +
  coord_cartesian(ylim = c(-2, 2)) +
  theme_bw() # 白色背景主题
dev.off()

library(pheatmap)
type = rep("HCC", ncol(res_combined))
names(type) = colnames(res_combined)
type[ colnames(data0)[grep("_NT", colnames(data0))] ] = "NT"
type[ colnames(data0)[grep("_DN", colnames(data0))] ] = "DN"
type[ colnames(data0)[grep("_T", colnames(data0))] ] = "eHCC"
type[ colnames(data2)[grep("\\.11A", colnames(data2))] ] = "NT"
annotation_col = data.frame(batch=as.factor(batch), type=as.factor(type) )
rownames(annotation_col) = colnames(res_combined)
pdf("Corrected_MCPcounter_heatmap.pdf", width=12, height=5)
pheatmap( pmin(pmax(t(scale(t(res_combined))), -2), 2), scale="none", annotation_col=annotation_col, show_colnames=F )
dev.off()

save.image("MCPcounter.Rdata")





