\name{gradient}
\alias{gradient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gradient of log likelihood }
\description{
'gradient' calculates the gradient of the log likelihood function.
}
\usage{
gradient(delta)
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
This function used for the 'optim' method. Some 'optim' methods need a
gradient and have to approximate it in every iteration step. This slows
down the process a lot. Therefore, an explicit gradient function increases
convergence speed and ensures correct convergence even if the log 
likelihood function is not continuous due to the indicator function of 
the ZI parameter.

The design matrices have to be defined as X (for mean), W (for overdispersion)
and Z (for ZI).

n has to be defined as the number of observations dim(X)[1].

k.beta has to be defined as the length of beta.

k.alpha has to be defined as the length of alpha.

k.gamma has to be defined as the length of gamma.

t.i has to be defined as the exposure.
}
\examples{
Y <- c(3,0,2)
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
gradient(delta)
#[1] 423.8616107 704.1026859   7.5116293  -0.9998638
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
