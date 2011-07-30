\name{beta.sample}
\alias{beta.sample}

\title{
Resampling multiple-site dissimilarity for n sites
}
\description{
Resamples the 3 multiple-site dissimilarities (beta.SIM, beta.NES and beta.SOR) 
for a subset of sites of the original data frame.
}
\usage{
beta.sample(x, sites=10, samples=100)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species}
\item{sites}{ number of sites for which multiple-site dissimilarities will be computed}
\item{samples}{ number of repetitions}
}


\value{
The function returns a list with the following components:
\item{sampled.values}{ dataframe containing beta.SIM, beta.NES and beta.SOR for all samples}
\item{mean.values}{ vector containing the mean values of beta.SIM, beta.NES and beta.SOR among samples}
\item{sd.values}{ vector containing the sd values of beta.SIM, beta.NES and beta.SOR among samples}
}

\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. 
Global Ecology and Biogeography 19:134-143
}

\author{
Andrés Baselga and David Orme
}


\seealso{
\code{\link{beta.multi}}, \code{\link{beta.sample}}, \code{\link{beta.temp}}
}
\examples{
# Read the data for Northern and Southern European cerambycids
data(ceram.s)
data(ceram.n)

# Resample 100 times the multiple-site dissimilarities
# for 10 countries.
beta.ceram.s<-beta.sample(ceram.s, 10, 100)
beta.ceram.n<-beta.sample(ceram.n, 10, 100)

# Plot the distributions of beta.SIM in Southern Europe (red) 
# and Northern Europe (blue)
hist(beta.ceram.s$sampled.values$beta.SIM, , breaks=10, col="red", xlim=c(0,1))
hist(beta.ceram.n$sampled.values$beta.SIM, , breaks=10, col="blue", add=TRUE)
}

