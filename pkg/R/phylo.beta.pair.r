phylo.beta.pair<-function (x, tree, index.family = "sorensen")
{
    index.family <- match.arg(index.family, c("jaccard", "sorensen"))
    pbc <- phylo.betapart.core(x,tree)

        ############ Paired matrix to distance matrix conversion (utility function) #######################

    dist.mat <- function(com,pair) {

      ncom <- nrow(com)
      distmat <- matrix(nrow=ncom,ncol=ncom,0,dimnames=list(rownames(com),rownames(com)))
      st <- c(0,cumsum(seq(ncom-1,2)))+1
      end <- cumsum(seq(ncom-1,1))
      for (i in 1:(ncom-1)) distmat[i,(ncom:(seq(1,ncom)[i]))]=c(pair[end[i]:st[i]],0)
      distmat <- as.dist(t(distmat))
      return(distmat)

    } # end of function dist.mat

    switch(index.family, sorensen = {
        phylo.beta.sim <- pbc$min.not.shared/(pbc$min.not.shared + pbc$shared)

        phylo.beta.sne <- ((pbc$max.not.shared - pbc$min.not.shared)/((2 * pbc$shared) + pbc$sum.not.shared)) * (pbc$shared/(pbc$min.not.shared + pbc$shared))

        phylo.beta.sor <- pbc$sum.not.shared/(2 * pbc$shared + pbc$sum.not.shared)

        phylo.pairwise <- list(phylo.beta.sim = dist.mat(x,phylo.beta.sim), phylo.beta.sne = dist.mat(x,phylo.beta.sne), phylo.beta.sor = dist.mat(x,phylo.beta.sor))
    								},

    					 jaccard = {
        phylo.beta.jtu <- (2 * pbc$min.not.shared)/((2 * pbc$min.not.shared) + pbc$shared)

        phylo.beta.jne <- ((pbc$max.not.shared - pbc$min.not.shared)/(pbc$shared + pbc$sum.not.shared)) * (pbc$shared/((2 * pbc$min.not.shared) + pbc$shared))

        phylo.beta.jac <- pbc$sum.not.shared/(pbc$shared + pbc$sum.not.shared)

        phylo.pairwise <- list(phylo.beta.jtu = dist.mat(x,phylo.beta.jtu), phylo.beta.jne = dist.mat(x,phylo.beta.jne), phylo.beta.jac = dist.mat(x,phylo.beta.jac))
    								}

    ) # end of switch

    return(phylo.pairwise)

} # end of function