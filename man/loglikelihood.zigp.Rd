\name{loglikelihood.zigp}
\alias{loglikelihood.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Log Likelihood of ZIGP distribution }
\description{
'loglikelihood.zigp' calculates the log likelihood function of the ZIGP 
distribution.
}
\usage{
loglikelihood.zigp(delta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{ a parameter vector of length (p+r+q) with dim(X) = (n x p), 
          dim(W) = (n x r), dim(Z) = (n x q).  Create delta by pasting 
          'delta <- c(beta, alpha, gamma)'.
          beta is the vector of regression parameters for the mean modelling.
          alpha is the vector of regression parameters for overdispersion modelling.
          gamma is the vector of regression parameters for the ZI modelling.
}
}
\details{
The response has to be defined as Ysave.
The design matrices have to be defined as Xsave (for mean), Wsave (for overdispersion)
and Zsave (for ZI).

n has to be defined as the number of observations.

k.beta has to be defined as the length of beta.

k.alpha has to be defined as the length of alpha.

k.gamma has to be defined as the length of gamma.

Y has to be defined as the response vector.

t.i has to be defined as the exposure.
}
\examples{
Ysave <- c(3,0,2)
Xsave <- matrix(c(1:3,4,3,5),3,2)
Wsave <- c(3,-4,-1)
Zsave <- rep(1,3)
n <- dim(Xsave)[1]
t.i <- rep(1,n)
beta  <- c(5,-2)
alpha <- 3.4
gamma <- -10
k.beta <- length(beta)
k.alpha <- 1
k.gamma <- length(gamma)
delta <- c(beta,alpha,gamma)
loglikelihood.zigp(delta)
#[1] 160.2377
}
\keyword{ distribution }% at least one, from doc/KEYWORDS

