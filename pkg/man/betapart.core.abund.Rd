\name{betapart.core.abund}
\alias{betapart.core.abund}
\encoding{utf8}

\title{
Core calculations of betapart for abundance-based dissimilarity measures
}
\description{
Computes the basic quantities needed for computing the abundance-based multiple-site 
dissimilarity measuresc and pairwise dissimilarity matrices.
}
\usage{
betapart.core.abund(x)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species} 
}

\value{
The function returns an object of class \code{betapart.abund} with the following elements:
\item{multiple.shared.abund}{ the multiple-site interesection component in terms of abundances (AM)}
\item{pair.shared.abund}{ a matrix containing the agreement in abundance between pairs of sites (A)}
\item{min.not.shared.abund}{ a matrix containing the minimum disagreement in abundance between pairs of sites: min(B,C)}
\item{max.not.shared.abund}{ a matrix containing the maximum disagreement in abundance between pairs of sites between pairs of sites: max(B,C)}
\item{pair.not.shared.abund}{ a matrix containing the total disagreement in abundance between pairs of sites: B+C}
}

\references{
Baselga, A. 2013. Separating the two components of abundance-based dissimilarity: balanced changes in 
abundance vs. abundance gradients. Methods in Ecology and Evolution 4: 552–557

Legendre, P. 2014. Interpreting the replacement and richness difference components of beta diversity. 
Global Ecology and Biogeography, 23: 1324–1334

Baselga, A. 2016. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution, in press
}

\author{
Andrés Baselga
}


\seealso{
\code{\link{beta.multi.abund}}, \code{\link{beta.pair.abund}}, \code{\link{beta.sample.abund}}, \code{\link{betapart.core}}
}
\examples{
data(ceram.s)
ceram.core.s<-betapart.core.abund(ceram.s)
ceram.dist.bray<-beta.pair.abund(ceram.core.s, index.family="bray")
ceram.dist.ruz<-beta.pair.abund(ceram.core.s, index.family="ruz")
ceram.multi.bray<-beta.multi.abund(ceram.core.s, index.family="bray")
ceram.multi.ruz<-beta.multi.abund(ceram.core.s, index.family="ruz")

}
