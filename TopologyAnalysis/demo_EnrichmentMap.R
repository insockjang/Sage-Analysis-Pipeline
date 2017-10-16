require(devtools)
source_url("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/TopologyAnalysis/draw_network_cytoscape_RCy3.R")

# input processing for enrichment map
# gsea results that are summarized from demo_gsea.R
gseaOutcome <- read.delim("/gpfs/archive/RED/isjang/Project_CC90009/SolidTumor/PCOR/WholeSamples/gsea_REACTOME.txt", stringsAsFactors = F)

#filter for FDR and No_genes (30<=x<=300)
cutoff.fdr <- 0.05
id <- which(gseaOutcome$fdr <= cutoff.fdr & gseaOutcome$No_genes >=30 & gseaOutcome$No_genes <= 300)
reference <- gseaOutcome[id,]

# Run Cytoscape (>3.0) apps in your background
# Once you already setup to link Cytoscape, you will get "vignette" network
# if you ran your gsea from MsigDB, you have to select mat = sim.combine.msigdb
# if you ran your gsea from Graphite, you have to select mat = sim.combine.graphite
CW <- drawEnrichMapCytoscape(mat = sim.combine.graphite,NODE = reference$geneset,size = reference$No_genes,stats = -log10(reference$fdr) * sign(reference$nes),similarity.cutoff = 0.3)

# After you manually modify your layout, you can directly save it into png format
saveImage(CW, "EnrichmentMap_wholeSamples", 'png', 6400)