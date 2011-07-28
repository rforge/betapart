\name{beta.pair}
\alias{beta.pair}
\alias{betapart}

\title{
Pair-wise dissimilarities
}
\description{
Computes 3 distance matrices accounting for the spatial turnover (beta.sim) 
and nestedness components (beta.nes) of beta diversity, and the sum of
both values (beta.sor)
}
\usage{
beta.pair(x)
}

\arguments{
\item{x}{data frame, where rows are sites and columns are species. Alternatively {\code{x} can be a {\code{betapart} object derived from the {\code{betapart.core} function}
}

\value{
The function returns a list with three dissimilarity matrices:
\item{beta.sim}{\code{dist} object, dissimilarity matrix accounting for spatial turnover, measured as Simpson pair-wise dissimilarity}
\item{beta.nes}{\code{dist} object, dissimilarity matrix accounting for nestedness-resultant dissimilarity, measured as nestedness-resultant pair-wise dissimilarity}
\item{beta.sim}{\code{dist} object, dissimilarity matrix accounting for beta diversity, measured as Sorensen pair-wise dissimilarity}
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143
}

\author{
Andrés Baselga and David Orme
}


\seealso{
\code{\link{beta.multi}}, \code{\link{beta.sample}}
}
\examples{
data(ceram.s)
ceram.dist<-beta.pair(ceram.s)
}

