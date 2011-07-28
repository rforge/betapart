\name{betapart.core}
\alias{betapart.core}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
betapart.core(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
    x <- as.matrix(x)
    shared <- x \%*\% t(x)
    not.shared <- abs(sweep(shared, 2, diag(shared)))
    sumSi <- sum(diag(shared))
    St <- sum(colSums(x) > 0)
    a <- sumSi - St
    sum.not.shared <- not.shared + t(not.shared)
    max.not.shared <- pmax(not.shared, t(not.shared))
    min.not.shared <- pmin(not.shared, t(not.shared))
    computations <- list(sumSi = sumSi, St = St, a = a, shared = shared, 
        not.shared = not.shared, sum.not.shared = sum.not.shared, 
        max.not.shared = max.not.shared, min.not.shared = min.not.shared)
    class(computations) <- "betapart"
    return(computations)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
