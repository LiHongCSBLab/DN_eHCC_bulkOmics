######################################### oncoprint ############################################
library(ComplexHeatmap)
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "white", col = "gray"))
    },
    "1" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FFA500", col = NA)),
    "2" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#687942", col = NA)),
    "3" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E2BF63", col = NA)),
    "4" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#1EAA39", col = NA)),
    "5" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E51373", col = NA)),
    "6" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA )),
    "7" = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "blue", col = NA ))
)
NinN_CFG = read.table("NinN_CFG.txt", header=T, sep="\t", row.names=1)
tmpP = gsub("_.*", "", colnames(NinN_CFG))
tmpType = gsub("\\d+", "", gsub(".*_", "", gsub("_R\\d+", "", colnames(NinN_CFG))))
sample_annotation = HeatmapAnnotation( patient=anno_block(labels=unique(tmpP), gp=gpar(fill="lightgray"), labels_gp=gpar(col="black") ), 
  type=gsub("\\d+", "", gsub(".*_", "", gsub("_R\\d+", "", colnames(NinN_CFG)))), 
  col = list( type=c("DN"="lightblue", "Ca"="pink") )
  )
tmpOnco = oncoPrint( NinN_CFG,  top_annotation=sample_annotation, alter_fun = alter_fun, column_split=factor(tmpP, levels=unique(tmpP) ), column_order=colnames(NinN_CFG) , show_column_names=T )
pdf("NinN_CFG.oncoPrint.pdf", width=10)
tmpOnco
dev.off()

tmpG = rownames(NinN_CFG)[row_order(draw(tmpOnco))]
tmpG = tmpG[length(tmpG):1]
tmpG_sampleN = cbind(DN=apply(NinN_CFG[tmpG, tmpType=="DN"], 1, function(x){ sum(!is.na(x)& x!="") }),
 HCC=apply(NinN_CFG[tmpG, tmpType=="Ca"], 1, function(x){ sum(!is.na(x) & x!="") }) )
tmpG_sampleN[, "DN"] = tmpG_sampleN[, "DN"]/19*100
tmpG_sampleN[, "HCC"] = tmpG_sampleN[, "HCC"]/23*100
pdf("NinN_CFG.barplot.pdf")
par(mar=c(3,5,3,3))
barplot( t(tmpG_sampleN), horiz=TRUE, beside=TRUE, las=1, col=c("blue", "red") )
dev.off()

