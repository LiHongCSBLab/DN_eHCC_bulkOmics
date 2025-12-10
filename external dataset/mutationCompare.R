library(readxl)

NinN_CFG = as.matrix(read.table("NinN_CFG.txt", header=T, sep="\t", row.names=1))
genes = rownames(NinN_CFG)
tmpType = gsub("\\d+", "", gsub(".*_", "", gsub("_R\\d+", "", colnames(NinN_CFG))))
tmpP = gsub("_.*", "", colnames(NinN_CFG))

# 按sample统计
data0 = cbind(DN=apply(NinN_CFG[, tmpType=="DN"], 1, function(x){y=as.numeric(x); sum(!is.na(y)& y!=6 & y!=7) })[genes], HCC=apply(NinN_CFG[, tmpType=="Ca"], 1, function(x){y=as.numeric(x); sum(!is.na(y)& y!=6 & y!=7) })[genes] )
data0["TERT",] = c(3, 5)
n0 = 23

# 按patient统计
data0 = cbind(DN=apply(NinN_CFG[, tmpType=="DN"], 1, function(x){y=as.numeric(x); length(unique(tmpP[tmpType=="DN"][!is.na(y)& y!=6 & y!=7])) })[genes], HCC=apply(NinN_CFG[, tmpType=="Ca"], 1, function(x){y=as.numeric(x); length(unique(tmpP[tmpType=="DN"][!is.na(y)& y!=6 & y!=7])) })[genes] )
data0["TERT",] = c(3, 4)
n0 = 16

count0 = data0[,"HCC"]
names(count0) = rownames(data0)
freq0 = count0/n0
cbind( freq0, count0) 

data1 = as.data.frame(read_excel("38355797/Mutations_20250815.xlsx"))
n1 = length(unique(data1[,1]))
count1 = sapply(genes, function(x){
  length(unique((data1[data1[,"Gene"]==x & (data1[,"Classification"]!="5'UTR" & data1[,"Classification"]!="3'UTR" & data1[,"Classification"]!="promoter"), 1])))
})
count1["TERT"] = length(unique((data1[data1[,"Gene"]=="TERT" & (data1[,"Classification"]=="promoter" ), 1])))
freq1 = count1/n1
p1 = sapply(genes, function(x){
  fisher.test(matrix(c(count1[x], n1-count1[x], count0[x], n0-count0[x]),nrow=2) )$p.value
})
q1 = p.adjust(p1, method="BH")
cbind( freq1, count1, p1, q1 ) 

data2 = read.csv("cBioPortal_TCGA_363HCC_MutatedGenes.txt", header=T, sep="\t", row.names=1)
n2 = 363
count2 = data2[genes, "Number_of_MutatedSamples"]
names(count2) = genes
count2["MDM4"] = 0
freq2 = count2/n2
p2 = sapply(genes, function(x){
  fisher.test(matrix(c(count2[x], n2-count2[x], count0[x], n0-count0[x]),nrow=2) )$p.value
})
q2 = p.adjust(p2, method="BH")
cbind( freq2, count2, p2, q2 ) 

data3 = read.csv("cBioPortal_24798001_231HCC_MutatedGenes.txt", header=T, sep="\t", row.names=1)
n3 = 231
count3 = data3[genes, "Number_of_MutatedSamples"]
names(count3) = genes
count3["MYC"] = 0
freq3 = count3/n3
p3 = sapply(genes, function(x){
  fisher.test(matrix(c(count3[x], n3-count3[x], count0[x], n0-count0[x]),nrow=2) )$p.value
})
q3 = p.adjust(p3, method="BH")
cbind( freq3, count3, p3, q3 ) 

library(pheatmap)
pdf("mutationCompare.pdf")
result = cbind(freq0, freq1, freq2, freq3)
colnames(result) = c("our", "CLCA", "TCGA", "PMID2479800")
pheatmap( result, cluster_cols=F, scale="none", cluster_rows=F )
dev.off()

save.image("mutationCompare.Rdata")


