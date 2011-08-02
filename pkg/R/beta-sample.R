beta.sample<-function(x, index.family="index.family", sites=nrow(x), samples=1){
	
	pb <- txtProgressBar(min = 0, max = samples, style = 3)

	if (pmatch(index.family, "sorensen", nomatch = 0)){	

	# Create a matrix to save the results.

	results.n<-as.data.frame(matrix(nrow=samples, ncol=3))
	colnames(results.n)<-c("beta.SIM", "beta.SNE", "beta.SOR")

	# Loop on the selected number of samples

		for(i in 1:samples){

  		# Take a sample of the dataset with the specified number of sites 
  
  		position<-as.vector(1:nrow(x))
  		sample.position<-sample(position, sites)
  		x.sample<-x[sample.position,1:ncol(x)]
  
  		# Compute the three indices for this sample and save in the results matrix
  
  		x.beta<-beta.multi(x.sample, index.family)	

  		results.n[i,1]<-x.beta$beta.SIM
  		results.n[i,2]<-x.beta$beta.SNE
  		results.n[i,3]<-x.beta$beta.SOR
		
		Sys.sleep(0.1)
   		# update progress bar
   		setTxtProgressBar(pb, i)         
		}

	close(pb)
	result<-list(sampled.values=results.n,mean.values=mean(results.n), sd.values=sd(results.n))
	return(result)
	}

	if (pmatch(index.family, "jaccard", nomatch = 0)){	

	# Create a matrix to save the results.

	results.n<-as.data.frame(matrix(nrow=samples, ncol=3))
	colnames(results.n)<-c("beta.JTU", "beta.JNE", "beta.JAC")

	# Loop on the selected number of samples

		for(i in 1:samples){

  		# Take a sample of the dataset with the specified number of sites 
  
  		position<-as.vector(1:nrow(x))
  		sample.position<-sample(position, sites)
  		x.sample<-x[sample.position,1:ncol(x)]
  
  		# Compute the three indices for this sample and save in the results matrix
  
  		x.beta<-beta.multi(x.sample, index.family)	

  		results.n[i,1]<-x.beta$beta.JTU
  		results.n[i,2]<-x.beta$beta.JNE
  		results.n[i,3]<-x.beta$beta.JAC

		Sys.sleep(0.1)
   		# update progress bar
   		setTxtProgressBar(pb, i)         
		}

	close(pb)
	result<-list(sampled.values=results.n,mean.values=mean(results.n), sd.values=sd(results.n))
	return(result)

	}

	else{
	stop("invalid dissimilarity index family")
	}


}
