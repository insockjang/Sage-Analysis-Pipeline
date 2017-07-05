source("PathwayAnalysis/pathwayAnalysis_customized.R")

# run test with multicore with FET (fisher exact test)
results.fet<-pathwayAnalysis(pathwayFileName = "your pathwayGeneSet",Reference = B,Test.method = "FET",cores = 2)

# run test with multicore with GSEA
results.gsea<-pathwayAnalysis(pathwayFileName = "your pathwayGeneSet",Reference = myReference,Test.method = "GSEA",cores = 8)

# combine all p-value from GSEA
Pval<-c()
for(k in 1:length(results)){
  Pval<-c(Pval,results.gsea[[k]]$gseaResult$p.value)
}
Pval <- p.adjust(Pval,"BH")
data.frame(geneset = names(results.gsea), pvalue = Pval)

names(Pval)<-names(results.gsea)

# Do you want to plot all Enrichment score? then do following:
results[["your wanted GeneSet Name"]]$gseaPlot(myReference,results[["your wanted GeneSet Name"]]$gseaResult$geneset)
