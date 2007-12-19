\name{optimized.run}
\alias{optimized.run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculation of valid initial values }
\description{
'optimized.run' calculates valid initial values for the 'optim' method.
}
\usage{
optimized.run(Y, X, W, Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ response vector of length n. }
  \item{X}{ design matrix of dim (n x p) for mean modelling. }
  \item{W}{ design matrix of dim (n x r) for overdispersion modelling. }
  \item{Z}{ design matrix of dim (n x q) for zero inflation modelling. }
}
\details{
t.i has to be defined as the exposure.
n has to be defined as the number of observations dim(X)[1].
}
\examples{
Y <- c(3,0,2)
X <- matrix(c(1:3,4,3,5),3,2)
W <- c(3,-4,-1)
Z <- rep(1,3)
n <- dim(X)[1]
t.i <<- rep(1,n)
optimized.run(Y,X,W,Z)
#[1]           X1         X2          W          Z 
#[2]   -0.6783368  0.4559966  0.1238029 -1.3862944 
}
\keyword{ robust }% at least one, from doc/KEYWORDS

