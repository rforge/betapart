\name{betapart}
\alias{betapart-package}
\alias{betapart}
\docType{package}
\title{Partitioning beta diversity into turnover and nestedness components}
\description{
\pkg{betapart} allows computing pair-wise dissimilarities (distance matrices)
and multiple-site dissimilarities separating the turnover and
nestedness components of beta diversity as described in Baselga (2010). 
}
\details{
The pairwise function \command{beta.pair} yields 3 distance matrices accounting 
for the spatial turnover (beta.sim) and the nestedness components (beta.nes)
of beta-diversity. The third distance matrix accounts for the sum of both 
componentes (beta-sor), i.e. beta diversity. The multiple site function \command{beta.multi} 
yields the spatial turnover (beta.SIM) and the nestedness components (beta.NES) 
of overall beta diversity, and the sum of both components (beta.SOR), i.e. beta 
diversity. The multiple-site values can be randomly sampled a specified number of 
times for a specified number of sites using the function \command{beta.sample}.
}
\author{
Andrés Baselga and David Orme

Maintainer: Andrés Baselga <andres.baselga@usc.es>
}
\references{
Baselga, A. 2010. Partitioning the turnover and nestedness components of beta diversity. 
Global Ecology and Biogeography 19:134-143
}

\examples{

}
