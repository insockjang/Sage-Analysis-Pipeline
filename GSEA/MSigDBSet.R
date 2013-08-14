MSigDBSet<-function(pathwayName = NULL){

  ################################################################## bs analysis
  # MSigDB from synapse
  mSigDB_annotations <- loadEntity("syn105363")
  mSigDB_symbolID <- loadEntity("syn105350")
  DB<-mSigDB_symbolID$objects$MsigDB_symbolID  
  
  
  if(is.element(pathwayName,"BIOCARTA")){
    allPathways <- mSigDB_annotations$objects$C2$BIOCARTA
  }
  if(is.element(pathwayName,"KEGG")){
    allPathways <- mSigDB_annotations$objects$C2$KEGG
  }
  if(is.element(pathwayName,"REACTOME")){
    allPathways <- mSigDB_annotations$objects$C2$REACTOME
  }
  if(is.element(pathwayName,"GO_BP")){
    allPathways <- mSigDB_annotations$objects$C5$GO_BP
  }
  if(is.element(pathwayName,"GO_CC")){
    allPathways <- mSigDB_annotations$objects$C5$GO_CC
  }
  if(is.element(pathwayName,"GO_MF")){
    allPathways <- mSigDB_annotations$objects$C5$GO_MF
  }
  
  
  # preparing for FET : make total set
  geneAllSetList <-DB$genesets[is.element(DB$geneset.names,allPathways)]
    
  return(geneAllSetList)
}
