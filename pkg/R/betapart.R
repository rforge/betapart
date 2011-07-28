## multi site similarity values

beta.multi <- function(x){

	if (inherits(x, "betapart")){
		
		 maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
             minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])

                        # indices
                        beta.sim <- minbibj / (minbibj + x$a)
                        beta.nes <- (x$a / (minbibj + x$a)) * ((maxbibj - minbibj) / ((2 * x$a) + maxbibj + minbibj))
                        beta.sor <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * x$a))

            		multi <- list(beta.sim=beta.sim, beta.nes=beta.nes,beta.sor=beta.sor)

        	return(multi)
	}

	else{
		x<-as.matrix(x)
                shared <- x %*% t(x)
                not.shared <-  abs(sweep(shared, 2, diag(shared)))


                sum.not.shared <- not.shared + t(not.shared)
                max.not.shared <- pmax(not.shared, t(not.shared))
                min.not.shared <- pmin(not.shared, t(not.shared))
        
                
                sumSi <- sum(diag(shared)) # species by site richness
                St <- sum(colSums(x) > 0)  # regional species richness
                a <- sumSi - St            # multi site shared species term


                maxbibj <- sum(max.not.shared[lower.tri(max.not.shared)])
                minbibj <- sum(min.not.shared[lower.tri(min.not.shared)])


                        # indices
                        beta.sim <- minbibj / (minbibj + a)
                        beta.nes <- (a / (minbibj + a)) * ((maxbibj - minbibj) / ((2 * a) + maxbibj + minbibj))
                        beta.sor <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * a))

            		multi <- list(beta.sim=beta.sim, beta.nes=beta.nes,beta.sor=beta.sor)
        	return(multi)
	} 
}

## pairwise site similarity values
beta.pair <- function(x){
	if (inherits(x, "betapart")){
		
 		beta.sim <- x$min.not.shared / (x$min.not.shared + x$shared)
            beta.nes <- ((x$max.not.shared - x$min.not.shared) / ((2 * x$shared) + x$sum.not.shared)) * (x$shared / (x$min.not.shared + x$shared))
            beta.sor <- x$sum.not.shared / (2 * x$shared + x$sum.not.shared)
                
            pairwise <- list(beta.sim=as.dist(beta.sim), beta.nes=as.dist(beta.nes),beta.sor=as.dist(beta.sor))
                
            return(pairwise)
	}
	else{

		x<-as.matrix(x)
                shared <- x %*% t(x)
                not.shared <-  abs(sweep(shared, 2, diag(shared)))


                sum.not.shared <- not.shared + t(not.shared)
                max.not.shared <- pmax(not.shared, t(not.shared))
                min.not.shared <- pmin(not.shared, t(not.shared))
                
                        beta.sim <- min.not.shared / (min.not.shared + shared)
                        beta.nes <- ((max.not.shared - min.not.shared) / ((2 * shared) + sum.not.shared)) * (shared / (min.not.shared + shared))
                        beta.sor <- sum.not.shared / (2 * shared + sum.not.shared)
                
                        pairwise <- list(beta.sim=as.dist(beta.sim), beta.nes=as.dist(beta.nes),beta.sor=as.dist(beta.sor))
                
                return(pairwise)
        
        } 
}