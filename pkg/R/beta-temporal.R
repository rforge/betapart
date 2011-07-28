beta.temp <- function(x, y){

x<-as.matrix(x)
y<-as.matrix(y)

ai<-apply(x&y, 1, sum)
bi<-apply(x&!y, 1, sum)
ci<-apply(!x&y, 1, sum)

beta.sor<- (bi+ci) / (2*ai+bi+ci)
beta.sim<- pmin(bi,ci)/(ai+pmin(bi,ci))
beta.nes<- beta.sor - beta.sim

result<-data.frame(beta.sim, beta.nes, beta.sor)
return(result)
}
