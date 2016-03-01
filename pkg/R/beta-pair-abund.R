
beta.pair.abund<-function (x, index.family = "bray") 
{
    index.family <- match.arg(index.family, c("bray", "ruzicka"))
    if (!inherits(x, "betapart.abund")) {
        x <- betapart.core.abund(x)
    }
    switch(index.family, bray = {
        beta.bray.bal <- x$pair.min.not.shared/(x$pair.min.not.shared + x$pair.shared)
        beta.bray.gra <- ((x$pair.max.not.shared - x$pair.min.not.shared)/((2 * 
            x$pair.shared) + x$pair.sum.not.shared)) * (x$pair.shared/(x$pair.min.not.shared + 
            x$pair.shared))
        beta.bray <- x$pair.sum.not.shared/(2 * x$pair.shared + x$pair.sum.not.shared)
        pairwise <- list(beta.bray.bal = as.dist(beta.bray.bal), beta.bray.gra = as.dist(beta.bray.gra), 
            beta.bray = as.dist(beta.bray))
    }, ruzicka = {
        beta.ruz.bal <- (2 * x$pair.min.not.shared)/((2 * x$pair.min.not.shared) + 
            x$pair.shared)
        beta.ruz.gra <- ((x$pair.max.not.shared - x$pair.min.not.shared)/(x$pair.shared + 
            x$pair.sum.not.shared)) * (x$pair.shared/((2 * x$pair.min.not.shared) + 
            x$pair.shared))
        beta.ruz <- x$pair.sum.not.shared/(x$pair.shared + x$pair.sum.not.shared)
        pairwise <- list(beta.ruz.bal = as.dist(beta.ruz.bal), beta.ruz.gra = as.dist(beta.ruz.gra), 
            beta.ruz = as.dist(beta.ruz))
    })
    return(pairwise)
}
