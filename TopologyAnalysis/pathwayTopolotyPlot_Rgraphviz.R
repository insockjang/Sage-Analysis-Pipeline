library("Rgraphviz")
#`Unfolded Protein Response (UPR)`
g<-pathwayGraph(GRAPHITE$REACTOME$Apoptosis)

reference <- read.delim("/gpfs/archive/RED/isjang/Project_CC90009/partialCorrelation_conditioning_CRBN_tumortype.txt",header = F, stringsAsFactors = F)

# reference <- read.delim("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/OVARY//partialCorrelation_conditioning_CRBN_tumortype.txt",header = F, stringsAsFactors = F)

cc <- max(abs(reference$V2))
color.gradient <- function(x, colors=c("blue","white","red"), colsteps=256) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(-cc,cc, length.out=colsteps)) ] )
}


aa<-color.gradient(reference$V2)

NODE <- nodes(g)
id <- match(NODE,reference$V1)
C<- reference$V2[id]
COLOR <- aa[id]
COLOR[is.na(id)]<-"white"
names(COLOR)<-NODE

plot(reference$V2,col = aa, pch = 19, cex = 2)
plot(C,col = COLOR, pch = 19, cex = 3)


attrs <- getDefaultAttrs()
attrs$node$label <- NODE
attrs$node$fontsize <- length(NODE)

attrs$node$height<-length(NODE)
attrs$node$width<-length(NODE)



attrs$edge$arrowsize <-0.5
attrs$edge$color <-"grey"
attrs$edge$minlen <- 0

nAttrs <- list()
nAttrs$fillcolor <- COLOR


# layout might be # “dot”, “neato”, “twopi”, “circo”, or “fdp” osage

png("~/compare_apoptosis.png", width = 1200, height = 1200)
set.seed(1)
plot(g, nodeAttrs = nAttrs,attrs=attrs, "dot")
legend("topright",legend = c("Sensitive","Resistant"),bty = "n", pch = 19, cex = 3, col = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)))
# legend.col(COLOR,c(-cc, cc))
dev.off()



######################################################### gsea
load("/gpfs/archive/RED/isjang/Project_CC90009/pathwayAnalysis_PCOR_CC90009_Gene_given_CRBN_tumortype/All_SolidTumor/GSEA.pcor.rdata")
myReference <- reference$V2
names(myReference) <- reference$V1
results.gsea.REACTOME[["Unfolded Protein Response (UPR)"]]$gseaPlot(myReference,results.gsea.REACTOME[["Unfolded Protein Response (UPR)"]]$gseaResult$geneset)

legend.col <- function(col, lev){
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}