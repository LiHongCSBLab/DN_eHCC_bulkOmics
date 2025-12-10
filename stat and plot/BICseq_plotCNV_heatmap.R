library(copynumber)
library(GenomicRanges)

## combine multiple BIC-Seq results files into a segment matrix
getSeg_BICSeq<-function(files, geneGR=NULL, outfile){
  result = vector()
  for(file in files){
    print(file)
    sample = sub("\\..*", "", sub(".*\\/", "", file))
    data = read.table(file, header=T)
    data[,"chrom"] = as.character( data[,"chrom"] )
    for(i in 1:nrow(data)){
        arm = ""
        centromere_start = hg19chromosomes[hg19chromosomes[,"chr"]==data[i,"chrom"] & hg19chromosomes[,"label"]=="centromere", "start"]
        centromere_end = hg19chromosomes[hg19chromosomes[,"chr"]==data[i,"chrom"] & hg19chromosomes[,"label"]=="centromere", "end"]
        if( data[i,"end"] <= centromere_start ){ arm="p" }
        if( data[i,"start"] >= centromere_end ){ arm="q" }
        if( data[i,"start"]<centromere_start & data[i,"end"]>centromere_start & data[i,"end"]<centromere_end ){ arm="p"; data[i,"end"]=centromere_start }
        if( data[i,"end"]>centromere_end & data[i,"start"]>centromere_start & data[i,"start"]<centromere_end ){ arm="q"; data[i,"start"]=centromere_end }
        if( data[i,"start"]>centromere_start & data[i,"end"]<centromere_end ){ next }
        if( data[i,"start"]<centromere_start & data[i,"end"]>centromere_end ){ 
          result = c(result, sample, data[i,"chrom"], "p", data[i,"start"], centromere_start, "", data[i,"log2.copyRatio"], data[i, "pvalue"])
          result = c(result, sample, data[i,"chrom"], "q", centromere_end, data[i,"end"], "", data[i,"log2.copyRatio"], data[i, "pvalue"])
        }else{
          result = c(result, sample, data[i,"chrom"], arm, data[i,"start"], data[i,"end"], data[i,"binNum"], data[i,"log2.copyRatio"], data[i, "pvalue"])
        }
    }
  }
  result = t(matrix(result, nrow=8))
  result[,2] = gsub("chr", "", result[,2])
  seg = data.frame( sampleID=as.character(result[,1]), chrom=as.integer(result[,2]), arm=as.character(result[,3]), start.pos=as.numeric(result[,4]), end.pos=as.numeric(result[,5]), n.probes=as.numeric(result[,6]), mean=as.numeric(result[,7]), p=as.numeric(result[,8]), stringsAsFactors=F )
  
  if( !is.null(geneGR) ){
    segToGenes = t(apply(seg, 1, function(x){
      segGR = makeGRangesFromDataFrame( data.frame(chr=paste("chr", as.numeric(x["chrom"]), sep=""), start=as.numeric(x["start.pos"]), end=as.numeric(x["end.pos"]) ))
      overlapGenes = unique(as.character( geneLocation[ as.matrix(findOverlaps(geneGR, segGR))[, "queryHits"], "gene"] ))
      c(x, paste(overlapGenes, collapse=",") )
    }))
    colnames(segToGenes)[ncol(segToGenes)] = "genes"
    write.table( segToGenes, outfile, sep="\t", quote=F, row.names=F, col.names=T )
  }
  return(seg)
}

########################################################3
## hg19
hg19chromosomes = read.table("/meta/hg19/hg19chromosomes.txt", header=T)
hg19chromosomes[,"chr"] = as.character( hg19chromosomes[,"chr"] )
hg19chromosomes[,"label"] = as.character( hg19chromosomes[,"label"] )

## gene
geneLocation = read.table("/meta/ensemble_GRCh37.85/Homo_sapiens.GRCh37.85.protein_coding")
colnames(geneLocation) = c("chr", "start", "end", "strand", "gene")
geneGR = makeGRangesFromDataFrame( geneLocation )

## BIC-Seq
sampleInfor = as.matrix(read.table("sampleInfor.txt", header=T, row.names=1))
samples = rownames(sampleInfor)
seg = getSeg_BICSeq( sapply(samples, function(s){ paste(s, ".BICseq.SCNA", sep="") }), geneGR, "BIC-Seq.SCNA.SegToGenes" )
seg[,4]=as.numeric(seg[,4])
seg[,5]=as.numeric(seg[,5])
seg[,7]=as.numeric(seg[,7])
seg[,"sampleID"] = sampleInfor[ seg[,"sampleID"] , "finalID"]

pdf("BIC-Seq.SCNA.pdf", height=8, width=12)
seg = seg[ seg[,5]-seg[,4]>1e6, ]
plotAberration( seg, thres.gain=0.2, mar=c(5, 10, 4, 2) )
dev.off()



