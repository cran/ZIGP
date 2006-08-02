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
          gamma is the vector of regression parameters for the ZI modelling.}
}
\details{
The design matrices have to be defined as X (for mean), W (for overdispersion)
and Z (for ZI).

n has to be defined as the number of observations dim(X)[1].

k.beta has to be defined as the length of beta.

k.alpha has to be defined as the length of alpha.

k.gamma has to be defined as the length of gamma.

Y has to be defined as the response vector.

t.i has to be defined as the exposure.
}
\examples{
Y <- c(3,0,2)
X <- matrix(c(1:3,4,3,5),3,2)
W <- c(3,-4,-1)
Z <- rep(1,3)
n <- dim(X)[1]
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

## The function is currently defined as
function(delta)

{

  eta.mu <- double(n)

  eta.phi <- double(n)

  eta.omega <- double(n)

  s1 <- double(1)

  s2 <- double(1)



    if(k.beta == 1) {

      eta.mu <- X * delta[1]

    }

    else {

      beta <- delta[1:k.beta]

      eta.mu <- X \%*\% beta

    }

    

    if(k.alpha == 1) {

      eta.phi <- W * delta[k.beta + 1]

    }

    else {

      alpha <- delta[(k.beta + 1) : (k.beta + k.alpha)]

      eta.phi <- W \%*\% alpha

    }

    

    if(k.gamma == 1) {

      eta.omega <- Z * delta[k.beta + k.alpha + 1]

    }

    else {

      gamma <- delta[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)]

      eta.omega <- Z \%*\% gamma

    }





    mu.i <- t.i*exp(eta.mu)

    b.i <- exp(eta.phi)

    phi.i <- 1 + b.i

    k.i <- exp(eta.omega)

    temp1 <- (mu.i + b.i*Y)



# terms of log likelihood with 1{Y(i) = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* (log(k.i+exp(-1/phi.i*mu.i)) -

          log(1+k.i)))

# terms of log likelihood with 1{Y(i) > 0} (excluding factorial)

    s2 <- sum(ifelse(Y>0,1,0)*(-log(1+k.i)+log(t.i)+eta.mu+(Y-1)*log(temp1)-

          Y*log(phi.i)-1/phi.i*temp1) - lgamma(Y+1))



  l <- s1 + s2

  return( - l)

  }
}
\keyword{ distribution }% at least one, from doc/KEYWORDS

