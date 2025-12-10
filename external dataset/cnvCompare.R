
tmp = as.matrix(read.csv("sampleInfor.txt", header=T, sep="\t"))
tmp = tmp[tmp[,"group"]=="HCC", ]
sample2patient = tmp[,"patient"]
names(sample2patient) = tmp[,"sample"]

data0 = read.csv("BIC-Seq.SCNA.SegToGenes", header=T, sep="\t") 
data0 = data0[apply(data0, 1, function(x){ length(intersect(x["sampleID"], names(sample2patient)))>0 }), ]
data0_AMP = as.matrix(data0[data0[,"mean"]>0.2 & data0[,"p"]<0.05, ])
data0_AMP = t(matrix(unlist( apply(data0_AMP, 1, function(x){ y=unlist(strsplit(x["genes"], split=",")); as.vector(rbind(rep(x["sampleID"], length(y)), y)) }) ), nrow=2))
colnames(data0_AMP) = c("sample", "gene")
data0_DEL = as.matrix(data0[data0[,"mean"]< -0.2 & data0[,"p"]<0.05, ])
data0_DEL = t(matrix(unlist( apply(data0_DEL, 1, function(x){ y=unlist(strsplit(x["genes"], split=",")); as.vector(rbind(rep(x["sampleID"], length(y)), y)) }) ), nrow=2))
colnames(data0_DEL) = c("sample", "gene")

# 按sample统计
n0 = length(unique(data0[,1]))
data0_AMPfreq = tapply(1:nrow(data0_AMP), data0_AMP[,"gene"], function(x){ length(unique(data0_AMP[x, "sample"])) })/n0
data0_DELfreq = tapply(1:nrow(data0_DEL), data0_DEL[,"gene"], function(x){ length(unique(data0_DEL[x, "sample"])) })/n0

# 按patient统计
n0 = length(unique(sample2patient[data0[,1]]))
data0_AMPfreq = tapply(1:nrow(data0_AMP), data0_AMP[,"gene"], function(x){ length(unique(sample2patient[data0_AMP[x, "sample"]])) })/n0
data0_DELfreq = tapply(1:nrow(data0_DEL), data0_DEL[,"gene"], function(x){ length(unique(sample2patient[data0_DEL[x, "sample"]])) })/n0


########################
library(readxl)
library(GenomicRanges)
geneLocation = read.table("Homo_sapiens.GRCh37.85.protein_coding")
colnames(geneLocation) = c("chr", "start", "end", "strand", "gene")
geneGR = makeGRangesFromDataFrame( geneLocation )
data1_seg = as.data.frame(read_excel("38355797/Copy_Number_Alteration_20250815.xlsx"))
data1_genes = t(apply(data1_seg, 1, function(x){
      segGR = makeGRangesFromDataFrame( data.frame(chr=x["Chr"], start=as.numeric(x["Start"]), end=as.numeric(x["End"]) ))
      overlapGenes = unique(as.character( geneLocation[ as.matrix(findOverlaps(geneGR, segGR))[, "queryHits"], "gene"] ))
      c(x, paste(overlapGenes, collapse=",") )
    }))
colnames(data1_genes)[ncol(data1_genes)] = "genes"
n1 = length(unique(data1_seg[,1]))
data1_AMP = as.matrix(data1_genes[as.numeric(data1_genes[,"CopyNumber"])>2^1.2, ])
data1_AMP = t(matrix(unlist( apply(data1_AMP, 1, function(x){ y=unlist(strsplit(x["genes"], split=",")); as.vector(rbind(rep(x["CaseID"], length(y)), y)) }) ), nrow=2))
colnames(data1_AMP) = c("sample", "gene")
data1_AMPfreq = tapply(1:nrow(data1_AMP), data1_AMP[,"gene"], function(x){ length(unique(data1_AMP[x, "sample"])) })/n1
data1_DEL = as.matrix(data1_genes[as.numeric(data1_genes[,"CopyNumber"])<2^(-1.2), ])
data1_DEL = t(matrix(unlist( apply(data1_DEL, 1, function(x){ y=unlist(strsplit(x["genes"], split=",")); as.vector(rbind(rep(x["CaseID"], length(y)), y)) }) ), nrow=2))
colnames(data1_DEL) = c("sample", "gene")
data1_DELfreq = tapply(1:nrow(data1_DEL), data1_DEL[,"gene"], function(x){ length(unique(data1_DEL[x, "sample"])) })/n1

########################
data2 = read.csv("cBioPortal_TCGA_364HCC_CNAgenes.txt", header=T, sep="\t")
n2 = 364
data2_AMPfreq = as.numeric(gsub("%","",data2[data2[,"CNA"]=="AMP","Freq"]))/100
names(data2_AMPfreq) = data2[data2[,"CNA"]=="AMP","Gene"]
data2_DELfreq = as.numeric(gsub("%","",data2[data2[,"CNA"]=="HOMDEL","Freq"]))/100
names(data2_DELfreq) = data2[data2[,"CNA"]=="HOMDEL","Gene"]

########################
data3 = read.csv("cBioPortal_24798001_231HCC_CNAgenes.txt", header=T, sep="\t")
n3 = 231
data3_AMPfreq = as.numeric(gsub("%","",data3[data3[,"CNA"]=="AMP","Freq"]))/100
names(data3_AMPfreq) = data3[data3[,"CNA"]=="AMP","Gene"]
data3_DELfreq = as.numeric(gsub("%","",data3[data3[,"CNA"]=="HOMDEL","Freq"]))/100
names(data3_DELfreq) = data3[data3[,"CNA"]=="HOMDEL","Gene"]

library(pheatmap)
pdf("cnvCompare.pdf")
genes = rownames(read.table("CFG.txt", header=T, sep="\t", row.names=1))
genes_AMP = intersect(intersect(intersect(names(data0_AMPfreq), names(data2_AMPfreq)), names(data3_AMPfreq)), genes)
tmp = cbind(our=data0_AMPfreq[genes_AMP], CLCA=data1_AMPfreq[genes_AMP], TCGA=data2_AMPfreq[genes_AMP], PMID2479800=data3_AMPfreq[genes_AMP])
pheatmap( tmp, cluster_cols=F, scale="none" )
write.csv(tmp, "CFG_AMP.csv")
genes_DEL = intersect(intersect(intersect(names(data0_DELfreq), names(data2_DELfreq)), names(data3_DELfreq)), genes)
tmp = cbind(our=data0_DELfreq[genes_DEL], CLCA=data1_DELfreq[genes_DEL], TCGA=data2_DELfreq[genes_DEL], PMID2479800=data3_AMPfreq[genes_DEL])
pheatmap( tmp, cluster_cols=F, scale="none" )
write.csv(tmp, "CFG_DEL.csv")
dev.off()

save.image("cnvCompare.Rdata")





