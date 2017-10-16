require(devtools)
require(graphite)
source_url("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/TopologyAnalysis/draw_network_cytoscape_RCy3.R")


# What you need is your interested pathway name and input reference for node specification
# e.g. "Unfolded Protein Respons(UPR)"
pathway.msigdb <- "REACTOME_UNFOLDED_PROTEIN_RESPONSE"

# Run Cytoscape (>3.0) apps in your background
# Once you already setup to link Cytoscape, you will get "vignette" network
# reference is the statistics two column data.frame: 1st column(V1) : gene symbol, 2nd column(V2) : gene statistics such as correlation value, fdr, etc
# artificially generate reference from MSigDB : REACTOME_UNFOLDED_PROTEIN_RESPONSE
reference <- runif(80)
names(reference)<-MSigDB$C2.CP.REACTOME$genesets[["REACTOME_UNFOLDED_PROTEIN_RESPONSE"]]

CW <- drawPathwayCytoscapeMSigDB(pathway.msigdb, reference)

# After you manually modify your layout, you can directly save it into png format
saveImage(CW, "UPR_demo_MSigDB from Graphite Topology", 'png', 6400)
