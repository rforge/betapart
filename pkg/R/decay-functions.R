########################################################################################
##### A function for fitting a distance-decay models

decay.model<-function(y, x, model.type="exponential", y.type="similarities", perm=100){
model.type <- match.arg(model.type, c("exponential", "power"))

switch(model.type, exponential = {

y.type <- match.arg(y.type, c("similarities", "dissimilarities"))

switch(y.type, similarities = {
y<-as.vector(y)
x<-as.vector(x)
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)

parameters <- list(data = data.frame(x,y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = exp(log.glm$coefficients[1]), b.slope = log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
return(c("Input data are similarities",parameters))

}, dissimilarities = {
y<-as.vector(1-y)
x<-as.vector(x)
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)


parameters <- list(data = data.frame(x,1-y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = 1-exp(log.glm$coefficients[1]), b.slope = -log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
return(c("Input data are dissimilarities",parameters))
})
}

, power = {

y.type <- match.arg(y.type, c("similarities", "dissimilarities"))


switch(y.type, similarities = {
y<-as.vector(y)
x<-as.vector(log(x))
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)

parameters <- list(data = data.frame(exp(x),y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = exp(log.glm$coefficients[1]), b.slope = log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
return(c("Input data are similarities",parameters))

}, dissimilarities = {
y<-as.vector(1-y)
x<-as.vector(log(x))
log.glm<-glm(y~x, family=gaussian(link = "log"))

null.dev<-vector(mode="numeric", length=perm)
for (i in 1:perm){
null.dev[i]<-glm(y~sample(x), family=gaussian(link = "log"))$deviance
}
p.value<-mean(null.dev<log.glm$deviance)


parameters <- list(data = data.frame(exp(x),1-y), model = log.glm, model.type = model.type, y.type = y.type, 
	pseudo.r.squared = 1-log.glm$deviance/log.glm$null.deviance, 
	a.intercept = 1-exp(log.glm$coefficients[1]), b.slope = -log.glm$coefficients[2], 
	p.value = ifelse(p.value==0,1/perm,p.value))
return(c("Input data are dissimilarities",parameters))

})
})
}

################################################################################




################################################################################
##### A function for plotting the decay models

plot.decay<-function(decay.model, xlim=c(0,max(decay.model$data[,1])), ylim=c(0,1), add=FALSE, remove.dots=FALSE, 
col="black", pch=1, lty=1, lwd=5, cex=1, ...){
if(!remove.dots){pch=pch}
else{pch=""}
dista<-sort(unique(decay.model$data[,1]))
model.type <- match.arg(decay.model$model.type, c("exponential", "power"))
y.type <- match.arg(decay.model$y.type, c("similarities", "dissimilarities"))

switch(model.type, exponential = {
if(!add){
switch(y.type, similarities = {
plot(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, decay.model$a.intercept*exp(decay.model$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
plot(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-decay.model$a.intercept)*exp(-decay.model$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
})
}
if(add){
switch(y.type, similarities = {
points(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, decay.model$a.intercept*exp(decay.model$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
points(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-decay.model$a.intercept)*exp(-decay.model$b.slope*dista), col=col, lty=lty, lwd=lwd, ...)
})
}
}
, power = {
if(!add){
switch(y.type, similarities = {
plot(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, decay.model$a.intercept*dista^decay.model$b.slope, col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
plot(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-decay.model$a.intercept)*dista^-decay.model$b.slope, col=col, lty=lty, lwd=lwd, ...)
})
}
if(add){
switch(y.type, similarities = {
points(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Similarity", col=col, pch=pch, cex=cex, ...)
lines(dista, decay.model$a.intercept*dista^decay.model$b.slope, col=col, lty=lty, lwd=lwd, ...)
}, dissimilarities = {
points(decay.model$data, xlim=xlim, ylim=ylim, xlab="Distance", ylab="Dissimilarity", col=col, pch=pch, cex=cex, ...)
lines(dista, 1-(1-decay.model$a.intercept)*dista^-decay.model$b.slope, col=col, lty=lty, lwd=lwd, ...)
})
}
})
}

################################################################################




################################################################################
##### A function to bootstrap the parameters
##### It takes a decay model (x) and bootstrap the parameters R times

boot.coefs.decay<-function(x, R){
ptm <- proc.time()
original.coefs<-c(x$a.intercept, x$b.slope)
names(original.coefs)<-c("a.intercept", "b.slope")
boot.coefs<-matrix(nrow=R, ncol=2)
colnames(boot.coefs)<-c("a.intercept", "b.slope")
for (i in 1:R){
data.sample<-x$data[sample(1:nrow(x$data), replace=TRUE), ]
boot.mod<-decay.model(data.sample[,2], data.sample[,1], model.type=x$model.type, y.type=x$y.type, perm=1)
boot.coefs[i,]<-c(boot.mod$a.intercept, boot.mod$b.slope)
}
mean.boot<-c(mean(boot.coefs[,1]), mean(boot.coefs[,2]))
names(mean.boot)<-c("a.intercept", "b.slope")
sd.boot=c(sd(boot.coefs[,1]), sd(boot.coefs[,2]))
names(sd.boot)<-c("a.intercept", "b.slope")

result<-list( model.type=x$model.type,y.type=x$y.type, original.coefs=original.coefs, boot.coefs=boot.coefs, mean.boot=mean.boot, sd.boot=sd.boot,time=proc.time() - ptm)
result
}

################################################################################
