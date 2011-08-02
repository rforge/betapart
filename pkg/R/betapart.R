## multi site similarity values

beta.multi <- function(x, index.family="index.family"){

	if (inherits(x, "betapart")){

	maxbibj <- sum(x$max.not.shared[lower.tri(x$max.not.shared)])
      minbibj <- sum(x$min.not.shared[lower.tri(x$min.not.shared)])
		
		if (pmatch(index.family, "sorensen", nomatch = 0)){	
		 
            # indices
            beta.sim <- minbibj / (minbibj + x$a)
            beta.sne <- (x$a / (minbibj + x$a)) * ((maxbibj - minbibj) / ((2 * x$a) + maxbibj + minbibj))
            beta.sor <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * x$a))

           	multi <- list(beta.SIM=beta.sim, beta.SNE=beta.sne,beta.SOR=beta.sor)

        	return(multi)
		}

		if (pmatch(index.family, "jaccard", nomatch = 0)){

		# indices
            beta.jtu <- (2*minbibj) / ((2*minbibj) + x$a)
            beta.jne <- (x$a / ((2*minbibj) + x$a)) * ((maxbibj - minbibj) / ((x$a) + maxbibj + minbibj))
            beta.jac <- (minbibj + maxbibj) / (minbibj + maxbibj + x$a)

           	multi <- list(beta.JTU=beta.jtu, beta.JNE=beta.jne, beta.JAC=beta.jac)

        	return(multi)
		}

		else{
		stop("invalid dissimilarity index.family")
		}
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

			if (pmatch(index.family, "sorensen", nomatch = 0)){

                  # indices
                  beta.sim <- minbibj / (minbibj + a)
                  beta.sne <- (a / (minbibj + a)) * ((maxbibj - minbibj) / ((2 * a) + maxbibj + minbibj))
                  beta.sor <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * a))

            	multi <- list(beta.SIM=beta.sim, beta.SNE=beta.sne, beta.SOR=beta.sor)
        		return(multi)
			}

			if (pmatch(index.family, "jaccard", nomatch = 0)){

		      # indices
                  beta.jtu <- (2*minbibj) / ((2*minbibj) + a)
                  beta.jne <- (a / ((2*minbibj) + a)) * ((maxbibj - minbibj) / ((a) + maxbibj + minbibj))
                  beta.jac <- (minbibj + maxbibj) / (minbibj + maxbibj + a)

            	multi <- list(beta.JTU=beta.jtu, beta.JNE=beta.jne, beta.JAC=beta.jac)

        		return(multi)
			}

			else{
			stop("invalid dissimilarity index family")
			}
	}
}

## pairwise site similarity values
beta.pair <- function(x, index.family="index.family"){
	if (inherits(x, "betapart")){
		
		if (pmatch(index.family, "sorensen", nomatch = 0)){
 		beta.sim <- x$min.not.shared / (x$min.not.shared + x$shared)
            beta.sne <- ((x$max.not.shared - x$min.not.shared) / ((2 * x$shared) + x$sum.not.shared)) * (x$shared / (x$min.not.shared + x$shared))
            beta.sor <- x$sum.not.shared / (2 * x$shared + x$sum.not.shared)
                
            pairwise <- list(beta.sim=as.dist(beta.sim), beta.sne=as.dist(beta.sne),beta.sor=as.dist(beta.sor))
                
            return(pairwise)
		}

		if (pmatch(index.family, "jaccard", nomatch = 0)){
 		beta.jtu <- (2*x$min.not.shared) / ((2*x$min.not.shared) + x$shared)
            beta.jne <- ((x$max.not.shared - x$min.not.shared) / (x$shared + x$sum.not.shared)) * (x$shared / ((2*x$min.not.shared) + x$shared))
            beta.jac <- x$sum.not.shared / (x$shared + x$sum.not.shared)
                
            pairwise <- list(beta.jtu=as.dist(beta.jtu), beta.jne=as.dist(beta.jne),beta.jac=as.dist(beta.jac))
                
            return(pairwise)
		}

		else{
		stop("invalid dissimilarity index family")
		}

	}
	else{
		x<-as.matrix(x)
                shared <- x %*% t(x)
                not.shared <-  abs(sweep(shared, 2, diag(shared)))


                sum.not.shared <- not.shared + t(not.shared)
                max.not.shared <- pmax(not.shared, t(not.shared))
                min.not.shared <- pmin(not.shared, t(not.shared))

			if (pmatch(index.family, "sorensen", nomatch = 0)){
                
                        beta.sim <- min.not.shared / (min.not.shared + shared)
                        beta.sne <- ((max.not.shared - min.not.shared) / ((2 * shared) + sum.not.shared)) * (shared / (min.not.shared + shared))
                        beta.sor <- sum.not.shared / (2 * shared + sum.not.shared)
                
                        pairwise <- list(beta.sim=as.dist(beta.sim), beta.sne=as.dist(beta.sne),beta.sor=as.dist(beta.sor))
                
                	return(pairwise)
        		}

			if (pmatch(index.family, "jaccard", nomatch = 0)){
 			beta.jtu <- (2*min.not.shared) / ((2*min.not.shared) + shared)
            	beta.jne <- ((max.not.shared - min.not.shared) / (shared + sum.not.shared)) * (shared / ((2*min.not.shared) + shared))
           		beta.jac <- sum.not.shared / (shared + sum.not.shared)
                
            	pairwise <- list(beta.jtu=as.dist(beta.jtu), beta.jne=as.dist(beta.jne),beta.jac=as.dist(beta.jac))
                
            	return(pairwise)
			}

			else{
			stop("invalid dissimilarity index family")
			}
        } 
}