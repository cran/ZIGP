\name{response.zigp}
\alias{response.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ZIGP response generator }
\description{
'response.zigp' generates a ZIGP response vector for a given setting
of linear predictors.
}
\usage{
response.zigp(X, W, Z, beta, alpha, gamma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ design matrix of dim (n x p) for modelling of mean }
  \item{W}{ design matrix of dim (n x r) for modelling of overdispersion }
  \item{Z}{ design matrix of dim (n x q)  for modelling of zero inflation }
  \item{beta}{ regression parameters for mean of length p }
  \item{alpha}{ regression parameters for overdispersion of length r }
  \item{gamma}{ regression parameters for zero inflation of length q }
}
\value{
Generates a ZIGP response vector of length n.
}
\examples{
X <- matrix(c(1:3,4,3,5),3,2)
W <- c(3,-4,-1)
Z <- rep(1,3)
beta  <- c(5,-2)
alpha <- 3.4
gamma <- -10
# set seed for random variable generator
set.seed(1)
response.zigp(X,W,Z,beta,alpha,gamma)
#[1]   0  55 132
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
