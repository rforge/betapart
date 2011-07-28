## multi site similarity values

beta.multi <- function(x){

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
                        beta.SIM <- minbibj / (minbibj + a)
                        beta.NES <- (a / (minbibj + a)) * ((maxbibj - minbibj) / ((2 * a) + maxbibj + minbibj))
                        beta.SOR <- (minbibj + maxbibj) / (minbibj + maxbibj + (2 * a))

            		multi <- list(beta.SIM=beta.SIM, beta.NES=beta.NES,beta.SOR=beta.SOR)
        	return(multi)
	} 

## pairwise site similarity values
beta.pair <- function(x){

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