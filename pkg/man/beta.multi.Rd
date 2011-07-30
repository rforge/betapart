\name{beta.multi}
\alias{beta.multi}

\title{
Multiple-site dissimilarities
}
\description{
Computes 3 multiple-site dissimilarities accounting for the overall spatial turnover (beta.SIM) 
and the overall nestedness components (beta.NES) of beta diversity, and the sum of
both values (beta.SOR)
}
\usage{
beta.multi(x)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species. Alternatively \code{x} can be a \code{betapart} object derived from the \code{betapart.core} function}
}
\value{
The function returns a list with the three multiple site dissimilarity values:
\item{beta.SIM}{ value of the overall turnover component, measured as Simpson dissimilarity}
\item{beta.NES}{ value of the overall nestedness component, measured as Nestedness-resultant dissimilarity}
\item{beta.SOR}{ value of the overall beta diversity measured, as Sorensen dissimilarity}
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. Global Ecology and Biogeography 19:134-143
}

\author{
Andrés Baselga and David Orme
}


\seealso{
\code{\link{beta.pair}}, \code{\link{beta.sample}}, \code{\link{beta.core}}, \code{\link{beta.temp}}
}
\examples{
data(ceram.s)
ceram.beta<-beta.multi(ceram.s)
}



