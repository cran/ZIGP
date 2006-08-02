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
Generates a ZIGP response vector of length n.}
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

## The function is currently defined as
function(X, W, Z, beta, alpha, gamma)

{

# get number of observations from X, W and Z

  if(is.matrix(X)) {

    nx <- dim(X)[1]

  }

  else {

    nx <- length(X)

  }

  if(is.matrix(W)) {

    nw <- dim(W)[1]

  }

  else {

    nw <- length(W)

  }

  if(is.matrix(Z)) {

    nz <- dim(Z)[1]

  }

  else {

    nz <- length(Z)

  }



  if (nx==nz&nx==nw) {

  n <- nx

  mu <- double(n)

  phi <- double(n)

  omega <- double(n)

  Y <- double(n)



  if(is.matrix(X)) {

    mu <- exp(X\%*\%beta)

  }

  else {

    mu <- exp(X*beta)

  }

  if(is.matrix(W)) {

    phi <- 1 + exp(W\%*\%alpha)

  }

  else {

    phi <- 1 + exp(W*alpha)

  }

  if(is.matrix(Z)) {

    omega <- exp(Z\%*\%gamma)/(1+exp(Z\%*\%gamma))

  }

  else {

    omega <- exp(Z*gamma)/(1+exp(Z*gamma))

  }

# create ZIGP response for all i

  for(i in 1:n) {

    Y[i] <- rzigp(1, mu[i], phi[i], omega[i])

  }

  return(Y)

  }

  else {

  print("X, W and Z do not have the same numbers of rows!")

  }

  }
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
