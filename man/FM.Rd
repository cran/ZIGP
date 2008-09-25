\name{FM}
\alias{FM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fisher Information }
\description{
'FM' calculates the (Expected) Fisher Information matrix.
}
\usage{
FM(beta, alpha, gamma, X, W, Z, Offset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ regression parameters for mean of length p }
  \item{alpha}{ regression parameters for overdispersion of length r }
  \item{gamma}{ regression parameters for zero inflation of length q }
  \item{X}{ design matrix of dim (n x p) for mean modelling. }
  \item{W}{ design matrix of dim (n x r) for overdispersion modelling. }
  \item{Z}{ design matrix of dim (n x q) for zero inflation modelling. }
  \item{Offset}{ exposure for individual observation lengths. Defaults to a vector of 1.
          The offset MUST NOT be in 'log' scale. }
}
\examples{
X <- matrix(c(1:3,4,3,5),3,2)
W <- c(3,-4,-1)
Z <- rep(1,3)
n <- 3
beta  <- c(5,-2)
alpha <- 3.4
gamma <- -10
FM(beta, alpha, gamma, X, W, Z)
#              [,1]          [,2]         [,3]          [,4]
#[1,]  1.469180e+03  2.412237e+03 6.249334e-03 -8.400641e-11
#[2,]  2.412237e+03  3.965799e+03 1.040260e-02 -3.360257e-10
#[3,]  6.249334e-03  1.040260e-02 2.101615e-03  2.520099e-10
#[4,] -8.400641e-11 -3.360257e-10 2.520099e-10  9.079162e-05
}
\keyword{ distribution }% at least one, from doc/KEYWORDS

