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
    
    print(names(MSigDB))
    
    allPathways <- MSigDB[[pathwayName]]$genesets
    
  }
  
  if(toupper(synID) == "CPDB"){
    
    githubURL <- "https://raw.githubusercontent.com/insockjang/Sage-Analysis-Pipeline/master/PathwayAnalysis/CPDB.Rdata"
    load(url(githubURL))
    
    print(names(CPDB))
    
    allPathways <- CPDB[[pathwayName]]
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
    
    print(names(GRAPHITE))
    
    allPathways <- curateDB(GRAPHITE[[pathwayName]]@entries)
    
  }
    
  # FET test
  if(is.element(Test.method,"FET")){
    allGenes <- unique(unlist(allPathways))
    
    testSet <- Reference
    
    AllGenes <- union(testSet,allGenes)
    
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
