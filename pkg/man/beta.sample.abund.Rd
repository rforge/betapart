\name{beta.sample}
\alias{beta.sample}
\encoding{utf8}

\title{
Resampling multiple-site dissimilarity for n sites
}
\description{
Resamples the 3 multiple-site dissimilarities (turnover, nestedness-resultant fraction, and overall 
beta diversity) for a subset of sites of the original data frame.
}
\usage{
beta.sample.abund(x, index.family="bray", sites=10, samples=100)
}

\arguments{
\item{x}{ data frame, where rows are sites and columns are species}
\item{index.family}{ family of dissimilarity indices, partial match of \code{"bray"} or \code{"ruzicka"}.}
\item{sites}{ number of sites for which multiple-site dissimilarities will be computed. If not specified, default is all sites.}
\item{samples}{ number of repetitions. If not specified, default is 1.}
}


\value{
The function returns a list with a dataframe with the resampled 3 multiple-site dissimilarities 
(balanced variation fraction, abundance-gradient fraction and overall dissimilarity; see \code{\link{beta.multi.abund}}), 
a vector with the respective means and a vector with the respective standard deviation.

For \code{index.family="bray"}:
\item{sampled.values}{ dataframe containing beta.BRAY.BAL, beta.BRAY.GRA and beta.BRAY for all samples}
\item{mean.values}{ vector containing the mean values of beta.BRAY.BAL, beta.BRAY.GRA and beta.BRAY among samples}
\item{sd.values}{ vector containing the sd values of beta.BRAY.BAL, beta.BRAY.GRA and beta.BRAY among samples}


For \code{index.family="ruzicka"}:
\item{sampled.values}{ dataframe containing beta.RUZ.BAL, beta.RUZ.GRA and beta.RUZ for all samples}
\item{mean.values}{ vector containing the mean values of beta.RUZ.BAL, beta.RUZ.GRA and beta.RUZ among samples}
\item{sd.values}{ vector containing the sd values of beta.RUZ.BAL, beta.RUZ.GRA and beta.RUZ among samples}
}

\references{
Baselga, A. 2016. Partitioning abundance-based multiple-site dissimilarity into components: balanced variation in abundance and abundance gradients. Methods in Ecology and Evolution, in press
}

\author{
Andrés Baselga
}


\seealso{
\code{\link{beta.multi.abund}}, \code{\link{beta.sample}}
}
\examples{
# Read the data for Northern and Southern European cerambycids
data(ceram.s)
data(ceram.n)

# Resample 100 times the multiple-site dissimilarities
# for 10 countries.
beta.ceram.s<-beta.sample(ceram.s, index.family="bray", sites=10, samples=100)
beta.ceram.n<-beta.sample(ceram.n, index.family="bray", sites=10, samples=100)

# Plot the distributions of beta.SIM in Southern Europe (red) 
# and Northern Europe (blue)
plot(density(beta.ceram.s$sampled.values$beta.SIM), col="red", xlim=c(0,1))
lines(density(beta.ceram.n$sampled.values$beta.SIM), col="blue")

# Compute the p-value of difference in beta.SIM between South and North 
# (i.e. the probability of finding in the North a higher value than 
# in the South)
p.value.beta.BRAY.BAl<-length(which(beta.ceram.s$sampled.values$beta.BRAY.BAL<
beta.ceram.n$sampled.values$beta.BRAY.BAL))/100

p.value.beta.SIM
# The result is 0 and we used 100 samples, so p<0.01
}

