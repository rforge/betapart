beta.sample<-function(x, sites=nrow(x), samples=1){

# Create a matrix to save the results. In this case we want to compute
# three indices (beta.SOR, beta.SIM and beta.NES) and n samples

results.n<-as.data.frame(matrix(nrow=samples, ncol=3))
colnames(results.n)<-c("beta.SIM", "beta.NES", "beta.SOR")

# Loop on the selected number of samples

for(i in 1:samples){

  # Take a sample of the dataset with the specified number of sites 
  
  position<-as.vector(1:nrow(x))
  sample.position<-sample(position, sites)
  x.sample<-x[sample.position,1:ncol(x)]
  
  # Compute the three indices for this sample and save the results in the
  # results matrix
  
  x.beta<-beta.multi(x.sample)	

  results.n[i,1]<-x.beta$beta.SIM
  results.n[i,2]<-x.beta$beta.NES
  results.n[i,3]<-x.beta$beta.SOR
}
result<-list(sampled.values=results.n,mean.values=mean(results.n), sd.values=sd(results.n))
return(result)
}
