require(devtools)
source_url("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/pathwayAnalysis.R")

# Reference (numeric vector with name)
load("/gpfs/archive/RED/isjang/Project_CC90009/SolidTumor/PCOR/pcor_wholeSolidTumor.Rdata")

# preparation : my test set for FET(fisher exact test; hyper geometric test)
myReference <- names(sort(pcor.crbn,decreasing = F)[1:100])

# if you installed already doParallel package, it will be okay to run with multicore running
results.fet.NCI <-pathwayAnalysis(synID="GRAPHITE",pathwayName = "NCI",Reference = myReference,Test.method = "FET",cores = 4)


# basic function to summarize all GSEA results into one table
# sorted by FDR first and then NES
summaryFET <- function(results.fet){
  Pval.fet<-c()
  for(k in 1:length(results.fet)){
    Pval.fet<-c(Pval.fet,results.fet[[k]]$fetResult$p.value)
  }
  FDR.fet <- p.adjust(Pval.fet,"BH")
  Result <- data.frame(geneset = as.character(names(results.fet)), pval = Pval.fet, fdr = FDR.fet)#[-idx.nan,]
  Result.final <- Result[order(Result[,3],Result[,2], decreasing = c(F,F) ),]
  Result.final
}

# all REACTOME analysis was summarized
tab.fet.NCI <-summaryFET(results.fet.NCI)

# save your table into tab delimited format for later use in enrichment map drawing
write.table(tab.fet.NCI, file = "fet_NCI.txt",row.names = F, col.names = T, quote = F, sep = "\t")

