reactomeAnalysis("~/newPredictiveModel/PRISM/ActArea/E/bs_Lasso//aggregated_metrics.txt", "~/newPredictiveModel/PRISM/ActArea/E/bs_Lasso/")
reactomeAnalysis("~/newPredictiveModel/CTRPv2/ActArea/E/bs_Lasso//aggregated_metrics.txt", "~/newPredictiveModel/CTRPv2/ActArea/E/bs_Lasso/")

reactomeAnalysis("~/newPredictiveModel/PRISM/ActArea/E/bs_ENet//aggregated_metrics.txt", "~/newPredictiveModel/PRISM/ActArea/E/bs_ENet/")
reactomeAnalysis("~/newPredictiveModel/CTRPv2/ActArea/E/bs_ENet//aggregated_metrics.txt", "~/newPredictiveModel/CTRPv2/ActArea/E/bs_ENet/")




load("/gpfs/archive/RED/isjang/Project_CC90009/partialCorrelation_conditioning_CRBN_tumortype.Rdata")
write.table(pcor.cmi, file  = "/gpfs/archive/RED/isjang/Project_CC90009/partialCorrelation_conditioning_CRBN_tumortype.txt", quote = F, sep = "\t", row.names = T, col.names = F)
reactomeAnalysis1("/gpfs/archive/RED/isjang/Project_CC90009/partialCorrelation_conditioning_CRBN_tumortype.txt", "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/ALL_SolidTumor/")



load("/gpfs/archive/RED/isjang/Project_CC90009/partialCorrelation_conditioning_CRBN_sitespecific.Rdata")

for(site in names(pcor.site)){
  a <- pcor.site[[site]]
  names(a)<-names(pcor.cmi)
  write.table(a, file  = paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/",site,"/partialCorrelation_conditioning_CRBN_tumortype.txt"), quote = F, sep = "\t", row.names = T, col.names = F)
  reactomeAnalysis1(paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/",site,"/partialCorrelation_conditioning_CRBN_tumortype.txt"), paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/",site,"/"))
}

for(site in names(pcor.site)){
  a <- cor.site[[site]]
  names(a)<-names(pcor.cmi)
  write.table(a, file  = paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/",site,"/correlation_tumortype.txt"), quote = F, sep = "\t", row.names = T, col.names = F)
  reactomeAnalysis1(paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/",site,"/correlation_tumortype.txt"), paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/",site,"/"))
}


cor.all <-read.csv("/gpfs/archive/RED/isjang/Project_CC90009/PRISM/Prism_correlation_table.csv")
write.table(cor.all[,c("X", "Pearson_Correlation")], file  = paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/All_SolidTumor/correlation_tumortype.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
reactomeAnalysis1(paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/All_SolidTumor/correlation_tumortype.txt"), paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/All_SolidTumor//"))


pathway.cor<-vector("list", length(pcor.site)+1)
names(pathway.cor) <- c("All_SolidTumor",names(pcor.site))
for(site in names(pcor.site)){
  pathway.cor[[site]] <- read.delim(paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/",site,"/reactomePathway.txt"), header = T, stringsAsFactors = F)
}
pathway.cor[["All_SolidTumor"]] <- read.delim(paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/All_SolidTumor//reactomePathway.txt"), header = T, stringsAsFactors = F)

reactome.name <- unique(unlist(lapply(pathway.cor, function(x){x$Description})))

NES <- c()
PVAL <- c()
FDR <- c()
for(site in names(pathway.cor)){
  id <- match(reactome.name, pathway.cor[[site]]$Description)
  kk1<-pathway.cor[[site]][id, c("NES"),drop = F]
  kk2<-pathway.cor[[site]][id, c("pvalue"),drop = F]
  kk3<-pathway.cor[[site]][id, c("p.adjust"),drop = F]
  names(kk1)<-names(kk2)<-names(kk3)<-site
  NES<- cbind(NES,as.matrix(kk1))
  PVAL<- cbind(PVAL,as.matrix(kk2))
  FDR<- cbind(FDR,as.matrix(kk3))
}

rownames(NES)<-rownames(PVAL)<-rownames(FDR)<-reactome.name
saveRDS(NES, file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/NES.rds")
saveRDS(PVAL, file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/PVAL.rds")
saveRDS(FDR, file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_COR_CC90009/FDR.rds")

######

pathway.pcor<-vector("list", length(pcor.site)+1)
names(pathway.pcor) <- c("All_SolidTumor",names(pcor.site))
for(site in names(pcor.site)){
  pathway.pcor[[site]] <- read.delim(paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/",site,"/reactomePathway.txt"), header = T, stringsAsFactors = F)
}
pathway.pcor[["All_SolidTumor"]] <- read.delim(paste0("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype//All_SolidTumor//reactomePathway.txt"), header = T, stringsAsFactors = F)

reactome.name <- unique(unlist(lapply(pathway.pcor, function(x){x$Description})))

NES <- c()
PVAL <- c()
FDR <- c()
for(site in names(pathway.pcor)){
  id <- match(reactome.name, pathway.cor[[site]]$Description)
  kk1<-pathway.pcor[[site]][id, c("NES"),drop = F]
  kk2<-pathway.pcor[[site]][id, c("pvalue"),drop = F]
  kk3<-pathway.pcor[[site]][id, c("p.adjust"),drop = F]
  names(kk1)<-names(kk2)<-names(kk3)<-site
  NES<- cbind(NES,as.matrix(kk1))
  PVAL<- cbind(PVAL,as.matrix(kk2))
  FDR<- cbind(FDR,as.matrix(kk3))
}

rownames(NES)<-rownames(PVAL)<-rownames(FDR)<-reactome.name

rownames(NES)<-rownames(PVAL)<-rownames(FDR)<-reactome.name
saveRDS(NES, file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/NES.rds")
saveRDS(PVAL, file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/PVAL.rds")
saveRDS(FDR, file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/FDR.rds")


NES <- readRDS(file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/NES.rds")
PVAL <- readRDS(file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/PVAL.rds")
FDR <- readRDS(file = "/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/FDR.rds")

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
heatmap.2(FDR.sign[id0,kk$colInd], col = hmcol, trace="none", margin=c(15, 35), cexRow = 1, cexCol = 0.75, key = F,Rowv = T, Colv = T, dendrogram = "row", keysize = 0.25)
