# Demo for GSEA
require(devtools)
source_url("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/pathwayAnalysis.R")

# Reference (numeric vector with name)
load("/gpfs/archive/RED/isjang/Project_CC90009/SolidTumor/PCOR/pcor_wholeSolidTumor.Rdata")
myReference <- pcor.crbn

# if you installed already doParallel package, it will be okay to run with multicore running
# Using graphite DB for reactome analysis
results.gsea.REACTOME <-pathwayAnalysis(synID="GRAPHITE",pathwayName = "REACTOME",Reference = myReference,Test.method = "GSEA",cores = 4)

# Using MSigDB v6 for Hallmark of Cancer DB analysis
results.gsea.HallmarkCancer <-pathwayAnalysis(synID="MSigDB",pathwayName = "H.ALL",Reference = myReference,Test.method = "GSEA",cores = 4)

# how to plot individual GSEA plot
# select your specified pathway name, e.g. "Unfolded Protein Response (UPR)"
# Don't forget you have to have "myReference" to draw GSEA plot!!!!
results.gsea.REACTOME[["Unfolded Protein Response (UPR)"]]$gseaPlot(myReference,results.gsea.REACTOME[["Unfolded Protein Response (UPR)"]]$gseaResult$geneset)


# basic function to summarize all GSEA results into one table
# sorted by FDR first and then NES
summaryGSEA <- function(results.gsea){
  Pval.gsea<-c()
  NES <-c()
  NoLE <-c()
  NoG <-c()
  for(k in 1:length(results.gsea)){
    Pval.gsea<-c(Pval.gsea,results.gsea[[k]]$gseaResult$p.value)
    NES<-c(NES,results.gsea[[k]]$gseaResult$nes)
    NoLE<-c(NoLE, length(results.gsea[[k]]$gseaResult$leadEdge))
    NoG<-c(NoG, length(results.gsea[[k]]$gseaResult$geneset))
  }
  Pval.gsea[Pval.gsea == 0]<- 1e-4
  FDR.gsea <- p.adjust(Pval.gsea,"BH")
  Result <- data.frame(geneset = as.character(names(results.gsea)), pval = Pval.gsea, fdr = FDR.gsea, nes = NES, No_leadingEdge = NoLE, No_genes = NoG)#[-idx.nan,]
  Result.final <- Result[order(Result[,4], Result[,3], decreasing = c(F,F) ),]
  Result.final
}

# all REACTOME analysis was summarized
tab.gsea.REACTOME <-summaryGSEA(results.gsea.REACTOME)

# save your table into tab delimited format for later use in enrichment map drawing
write.table(tab.gsea.REACTOME, file = "gsea_REACTOME.txt",row.names = F, col.names = T, quote = F, sep = "\t")
