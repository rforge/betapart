\name{beta.temp}
\alias{beta.temp}

\title{
Temporal change in community composition
}
\description{
Computes the dissimilarity for each locality between time 1 and time 2, considering 
the turnover (beta.sim) and nestedness components (beta.nes) of temporal change, 
and the sum of both values (beta.sor)
}
\usage{
beta.temp(x, y)
}

\arguments{
\item{x}{ data frame for time 1, where rows are sites and columns are species}
\item{y}{ data frame for time 2, where rows are sites and columns are species. 
\code{x} and \code{y} must contain exactly the same sites and species}
}

\value{
The function returns a data frame where rows are sites and columns are pairwise dissimilarity 
values between cell composition in time 1 and time 2: beta.sim, beta.nes, and beta.sor
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143
}

\author{
Andrés Baselga and David Orme
}


\seealso{
\code{\link{beta.multi}}, \code{\link{beta.pair}}, \code{\link{beta.sample}}, \code{\link{beta.core}}, 
}
\examples{
data(ceram.s)
ceram.t<-beta.temp(ceram.s, ceram.s)
}