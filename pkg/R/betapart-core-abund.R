betapart.core.abund <- function(x){
	
	# test for a numeric matrix or data.frame
	if(! is.matrix(x)){
		x<-as.matrix(x)
  	}
	
	if(! is.numeric(x))
    	stop("The data in x is not numeric.",call.=TRUE)
	
	# simple test for positive data
	xvals <-  unique(as.vector(x))
	if (any(xvals<0)) 
        stop("The table contains negative values: data should be abundances.", call. = TRUE)

	pair.shared<-matrix(nrow=nrow(x),ncol=nrow(x))
	rownames(pair.shared)<-rownames(x)
	colnames(pair.shared)<-rownames(x)

	pair.not.shared<-matrix(nrow=nrow(x),ncol=nrow(x))
	rownames(pair.not.shared)<-rownames(x)
	colnames(pair.not.shared)<-rownames(x)

	for(i in 1:nrow(x)) {
    	for(j in i:nrow(x)) {
		pair.shared[j,i]<-sum(pmin(x[i,],x[j,]))
		pair.not.shared[i,j]<-sum(x[i,])-sum(pmin(x[i,],x[j,]))
		pair.not.shared[j,i]<-sum(x[j,])-sum(pmin(x[i,],x[j,]))
		}
		}

	pair.shared<-as.dist(pair.shared)
	pair.max.not.shared<-pmax(as.dist(t(upper.tri(pair.not.shared)*pair.not.shared)), as.dist(pair.not.shared))
	pair.min.not.shared<-pmin(as.dist(t(upper.tri(pair.not.shared)*pair.not.shared)), as.dist(pair.not.shared))

	multiple.shared<-sum(x)-sum(apply(x,2,max))



	computations<-list(data=x, multiple.shared.abund=multiple.shared, pair.shared.abund=pair.shared, pair.min.not.shared.abund=pair.min.not.shared, 
		pair.max.not.shared.abund=pair.max.not.shared, pair.not.shared.abund=pair.min.not.shared+pair.max.not.shared)
    class(computations)<-"betapart.abund"

	return(computations)
} 
