betapart.core <- function(x){
	x<-as.matrix(x)
  
	if(any(x>1))
    	stop("The table contains values >1: data should be presence/absence.",call.=TRUE);

	shared <- x %*% t(x)
      not.shared <-  abs(sweep(shared, 2, diag(shared)))
		
	sumSi <- sum(diag(shared)) # species by site richness
      St <- sum(colSums(x) > 0)  # regional species richness
      a <- sumSi - St            # multi site shared species term


      sum.not.shared <- not.shared + t(not.shared)
      max.not.shared <- pmax(not.shared, t(not.shared))
      min.not.shared <- pmin(not.shared, t(not.shared))

 	computations<-list(sumSi=sumSi, St=St, a=a, shared=shared, not.shared=not.shared, sum.not.shared=sum.not.shared, max.not.shared=max.not.shared, min.not.shared=min.not.shared)
      class(computations)<-"betapart"  
	return(computations)
} 