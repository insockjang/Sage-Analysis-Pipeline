reactomeAnalysis <- function(inputReferenceFile, outputFolder){
  library(org.Hs.eg.db)
  library(DOSE)
  library(ReactomePA)
  library(clusterProfiler)
  library(ReactomePA)
  
  reference <- read.delim(inputReferenceFile,header = F, stringsAsFactors = F)
  
  gL <- reference$V2 - min(reference$V2, na.rm = T)
  names(gL)<-reference$V1
  geneList <- gL
  
  X <- sort(geneList, decreasing = T, index.return = T)
  
  de <- names(X$x)
  head(de)
  DE <- bitr(de, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  
  ## ------------------------------------------------------------------------
  
  id <- match(DE$SYMBOL, de)
  Y<-X$x[id]
  names(Y)<-DE$ENTREZID
  
  y <- gsePathway(Y, nPerm=1000,
                  minGSSize=30, maxGSSize = 300, pvalueCutoff = 1,
                  pAdjustMethod="BH", verbose=FALSE)
  
  res <- as.data.frame(y)
  head(res)
  write.table(res, file = paste0(outputFolder,"/reactomePathway.txt"),quote = F, sep = "\t", col.names = T, row.names = F)
  ## ----fig.height=16, fig.width=16, eval=FALSE-----------------------------
  png(paste0(outputFolder, "enrichmentMap.png"),,width = 800, height = 800)
  enrichMap(y,layout=igraph::layout.auto,vertex.label.cex = 0.75,fig.height = 24, fig.width = 24)
  dev.off()
  
  # ## ----fig.height=7, fig.width=10------------------------------------------
  # gseaplot(y, geneSetID = "R-HSA-69242")
  # 
  # ## ----fig.height=16, fig.width=16, eval=FALSE-----------------------------
  # viewPathway("G alpha (i) signalling events", readable=TRUE, foldChange=geneList, vertex.label.cex = 0.5)
  # 
  
}