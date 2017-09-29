# work flow or Demo using MSigDB or Graphite DB
require(graphite)
require(doParallel)
require(devtools)
source_url("https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/myPathwayAnalysis.R")

pathwayAnalysis<-function(synID=NULL,pathwayName = NULL,Reference = NULL,Test.method = c("FET","GSEA"),cores = 1){
  
  if(is.null(synID)==1){
    error("Please select which database you would like to use in your analysis: MSigDB(syn1681370) or Graphite(syn2135025)")
    break;
  }
  if(is.null(pathwayName)==1){
    error("Please select which pathway you would like to use in your analysis: KEGG, BioCarta, ...")
    break;
  }
  if(is.null(Reference)==1){
    error("Please provide your reference list to be tested in pathway analysis: list of genes or two column reference list (gene name and its statistics")
    break;
  }
  
  pathwayName = toupper(pathwayName)
  
  if(toupper(synID) == "MSIGDB"){
    
    githubURL <- "https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/MSigDBv6.Rdata"
    load(url(githubURL))
    
    if(is.element(pathwayName,"BIOCARTA")){
      allPathways <- MSigDB$C2.CP.BIOCARTA$genesets
    }
    if(is.element(pathwayName,"KEGG")){
      allPathways <- MSigDB$C2.CP.KEGG$genesets
    }
    if(is.element(pathwayName,"REACTOME")){
      allPathways <- MSigDB$C2.CP.REACTOME$genesets
    }
    if(is.element(pathwayName,"GO_BP")){
      allPathways <- MSigDB$C5.BP$genesets
    }
    if(is.element(pathwayName,"GO_CC")){
      allPathways <- MSigDB$C5.CC$genesets
    }
    if(is.element(pathwayName,"GO_MF")){
      allPathways <- MSigDB$C5.MF$genesets
    }
    if(is.element(pathwayName,"HALLMARK")){
      allPathways <- MSigDB$H.ALL$genesets
    }
  }
  if(toupper(synID) == "GRAPHITE"){
    
    curateDB <- function(DB){
      db <- vector("list",length(DB))
      for(k in 1:length(DB)){
        db[[k]]<-nodes(DB[[k]])
      }
      names(db) <- names(DB)
      return(db)
    }
    
    
    githubURL <- "https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/Graphite.Rdata"
    load(url(githubURL))
    
    if(is.element(pathwayName,"BIOCARTA")){
      allPathways <- curateDB(GRAPHITE$BIOCARTA@entries)
    }
    if(is.element(pathwayName,"KEGG")){
      allPathways <- curateDB(GRAPHITE$KEGG@entries)
    }
    if(is.element(pathwayName,"REACTOME")){
      allPathways <- curateDB(GRAPHITE$REACTOME@entries)
    }
    if(is.element(pathwayName,"NCI")){
      allPathways <- curateDB(GRAPHITE$NCI@entries)
    }
    if(is.element(pathwayName,"HUMANCYC")){
      allPathways <- curateDB(GRAPHITE$HUMANCYC@entries)
    }
    if(is.element(pathwayName,"PANTHER")){
      allPathways <- curateDB(GRAPHITE$PANTHER@entries)
    }
  }
    
  # FET test
  if(is.element(Test.method,"FET")){
    allGenes <-c()
    for (i in 1:length(allPathways)){
      allGenes<-union(allGenes,allPathways[[i]])
    }      
    testSet = Reference
    AllGenes<-union(testSet,allGenes)
    
    pathwayTest <-function(x){    
      curPathwayGenes <- allPathways[[x]]        
      pathwayAnalysis<-myPathwayAnalysis$new()
      pathwayAnalysis$fet(AllGenes,curPathwayGenes,testSet)
      return(pathwayAnalysis)      
    }    
    results<-mclapply(1:length(allPathways), function(x)pathwayTest(x),mc.cores= cores)     
    
  }else{ # GSEA test
    referenceSet<-sort(Reference, decreasing =TRUE, index.return =TRUE)
    pathwayTest <-function(x){    
      curPathwayGenes <- allPathways[[x]]        
      pathwayAnalysis<-myPathwayAnalysis$new()
      pathwayAnalysis$gsea(referenceSet$x,curPathwayGenes,np=10000,w =1)
      return(pathwayAnalysis)      
    }     
    results<-mclapply(1:length(allPathways), function(x)pathwayTest(x),mc.cores= cores)               
  }
  
  names(results)<-names(allPathways)  
  return(results)  
}
