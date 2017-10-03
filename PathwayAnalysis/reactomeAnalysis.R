  NES <- readRDS(file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/NES.rds")
PVAL <- readRDS(file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/PVAL.rds")
FDR <- readRDS(file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/FDR.rds")

id0 <- which(apply(FDR,1,function(x){length(which(x<0.05))}) >2)
id1 <- which(apply(NES,1,function(x){length(which(abs(x)>2))}) !=0)

id <- intersect(id1,id0)

FDR.sign <- sign(NES) * -log10(FDR)

library("RColorBrewer")
library("gplots")

hmcol = colorpanel(100,"red","black","green")
kk<-heatmap.2(FDR.sign[id,])
heatmap.2(FDR.sign[id,kk$colInd], col = hmcol, trace="none", margin=c(15, 35), cexRow = 1, cexCol = 0.75, key = F,Rowv = T, Colv = T, dendrogram = "row", keysize = 0.25)



id0 <- which(apply(FDR,1,function(x){length(which(x<0.05))}) >2)
kk<-heatmap.2(FDR.sign[id0,])
png("~/cor_pathway_heatmat.png",width = 1200, height = 2400)
heatmap.2(FDR.sign[id0,kk$colInd], col = hmcol, trace="none", margin=c(15, 35), cexRow = 1, cexCol = 0.75, key = F,Rowv = T, Colv = T, dendrogram = "row", keysize = 0.25)
dev.off()
