beta.temp <- function(x, y, index.family="index.family"){

	x<-as.matrix(x)
	y<-as.matrix(y)

	ai<-apply(x&y, 1, sum)
	bi<-apply(x&!y, 1, sum)
	ci<-apply(!x&y, 1, sum)

	if (pmatch(index.family, "sorensen", nomatch = 0)){
		beta.sor<- (bi+ci) / (2*ai+bi+ci)
		beta.sim<- pmin(bi,ci)/(ai+pmin(bi,ci))
		beta.sne<- beta.sor - beta.sim

	result<-data.frame(beta.sim, beta.sne, beta.sor)
	return(result)
	}

	if (pmatch(index.family, "jaccard", nomatch = 0)){
		beta.jac<- (bi+ci) / (ai+bi+ci)
		beta.jtu<- 2*pmin(bi,ci)/(ai+(2*pmin(bi,ci)))
		beta.jne<- beta.jac - beta.jtu

	result<-data.frame(beta.jtu, beta.jne, beta.jac)
	return(result)
	}

	else{
	stop("invalid dissimilarity index family")
	}
}
