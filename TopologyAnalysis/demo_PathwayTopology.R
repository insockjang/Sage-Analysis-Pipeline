require(devtools)
source_url("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/TopologyAnalysis/draw_network_cytoscape_RCy3.R")
githubURL <- "https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/Graphite.Rdata"
load(url(githubURL))


# What you need is your interested pathway name and input reference for node specification
# e.g. "Unfolded Protein Respons(UPR)"
pathway <- pathwayGraph(GRAPHITE$REACTOME[["Unfolded Protein Response (UPR)"]])

# Run Cytoscape (>3.0) apps in your background
# Once you already setup to link Cytoscape, you will get "vignette" network
# reference is the statistics two column data.frame: 1st column(V1) : gene symbol, 2nd column(V2) : gene statistics such as correlation value, fdr, etc
load("/gpfs/archive/RED/isjang/Project_CC90009/SolidTumor/PCOR/pcor_wholeSolidTumor.Rdata")

myReference <- data.frame(V1 = names(pcor.crbn), V2 = pcor.crbn)

CW <- drawPathwayCytoscapeAdvance(pathway, myReference)

# After you manually modify your layout, you can directly save it into png format
saveImage(CW, "UPR_pcor", 'png', 6400)
