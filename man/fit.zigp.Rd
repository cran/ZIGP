\name{fit.zigp}
\alias{fit.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitted Values }
\description{
'fit.zigp' computes fitted values.
}
\usage{
fit.zigp(delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{ a parameter vector of length (p+r+q) with dim(X) = (n x p), 
          dim(W) = (n x r), dim(Z) = (n x q).  Create delta by pasting 
          'delta <- c(beta, alpha, gamma)'.
          beta is the vector of regression parameters for the mean modelling.
          alpha is the vector of regression parameters for overdispersion modelling.
          gamma is the vector of regression parameters for the ZI modelling. }
}
\details{
The design matrices have to be defined as X (for mean), W (for overdispersion)
and Z (for ZI).

n has to be defined as the number of observations dim(X)[1].

k.beta has to be defined as the length of beta.

k.alpha has to be defined as the length of alpha.

k.gamma has to be defined as the length of gamma.

t.i has to be defined as the exposure.
}
\examples{
X <- matrix(c(1:3,4,3,5),3,2)
W <- c(3,-4,-1)
Z <- rep(1,3)
n <- dim(X)[1]
beta  <- c(5,-2)
alpha <- 3.4
gamma <- -10
k.beta <- length(beta)
k.alpha <- 1
k.gamma <- length(gamma)
t.i <- rep(1,n)
delta <- c(beta,alpha,gamma)
fit.zigp(delta)
#[1]   $fit
#[2]                [,1]
#[3]   [1,]   0.04978481
#[4]   [2,]  54.59567139
#[5]   [3,] 148.40642146
#[6]   $mu
#[7]                [,1]
#[8]   [1,]   0.04978707
#[9]   [2,]  54.59815003
#[10]  [3,] 148.41315910
#[11]  $phi
#[12]  [1] 26904.186074     1.000001     1.033373
#[13]  $omega
#[14]  [1] 4.539787e-05 4.539787e-05 4.539787e-05
}
\keyword{ datagen }% at least one, from doc/KEYWORDS

