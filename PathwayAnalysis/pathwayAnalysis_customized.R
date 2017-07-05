# work flow or Demo using MSigDB or Graphite DB
require(synapseClient)
require(parallel)
source("/home/isjang/Code/Sage-Analysis-Pipeline/PathwayAnalysis/myPathwayAnalysis.R")

pathwayAnalysis<-function(pathwayFileName = NULL,Reference = NULL,Test.method = c("FET","GSEA"),cores = 1){
  
  myPathwayGeneSet <- readRDS(pathwayFileName)
  allPathways <- myPathwayGeneSet
  
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
      pathwayAnalysis$gsea(referenceSet$x,curPathwayGenes,np=1000,w =1)
      return(pathwayAnalysis)      
    }     
    results<-mclapply(1:length(allPathways), function(x)pathwayTest(x),mc.cores= cores)               
  }
  
  names(results)<-names(allPathways)  
  return(results)  
}
