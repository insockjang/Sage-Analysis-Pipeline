# weighted aggregation
weightAggreation<-function(resultsModel){
  
  ResultBS<-c()
  for(k in 1:length(resultsModel)){  
    ResultBS<-cbind(ResultBS,rank(abs(as.numeric(resultsModel[[k]][-1])),ties.method="min")/length(resultsModel[[k]]))
  }
  rownames(ResultBS)<-rownames(resultsModel[[1]])[-1]  
  reference <- apply(ResultBS,1,sum) 
  reference<-sort(reference,decreasing = T)
  return(reference)
}

# linear(count) aggregation
countAggreation<-function(resultsModel){
  
  ResultBS<-c()
  for(k in 1:length(resultsModel)){    
    ResultBS<-cbind(ResultBS,as.matrix(resultsModel[[k]]!=0))
  }
  ResultBS<-ResultBS[-1,]
  rownames(ResultBS)<-rownames(resultsModel[[1]])[-1]  
  reference <- apply(ResultBS,1,sum) 
  reference<-sort(reference,decreasing = T)
  
  return(reference)
}

averageAggreation<-function(resultsModel){  
  ResultBS<-c()
  for(k in 1:length(resultsModel)){  
    ResultBS<-cbind(ResultBS,resultsModel[[k]][-1])
  }
  rownames(ResultBS)<-rownames(resultsModel[[1]])[-1]  
  reference <- abs(apply(ResultBS,1,mean))
  
  reference<-sort(reference,decreasing = T)
  return(reference)
}

variationAggreation<-function(resultsModel){  
  ResultBS<-c()
  for(k in 1:length(resultsModel)){  
    ResultBS<-cbind(ResultBS,resultsModel[[k]][-1])
  }
  rownames(ResultBS)<-rownames(resultsModel[[1]])[-1]  
  reference <- apply(ResultBS,1,sd) 
  reference<-sort(reference,decreasing = T)
  return(reference)
}


absoluteAverageAggreation<-function(resultsModel){  
  ResultBS<-c()
  for(k in 1:length(resultsModel)){  
    ResultBS<-cbind(ResultBS,abs(resultsModel[[k]][-1]))
  }
  rownames(ResultBS)<-rownames(resultsModel[[1]])[-1]  
  reference <- abs(apply(ResultBS,1,mean))
  
  reference<-sort(reference,decreasing = T)
  return(reference)
}
