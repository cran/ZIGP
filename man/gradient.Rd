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

## The function is currently defined as
function(delta)

{

  grad <- double(k.beta + k.alpha + k.gamma)

  eta.mu <- double(n)

  eta.phi <- double(n)

  eta.omega <- double(n)

  s1 <- double(1)

  s2 <- double(1)



    if(k.beta == 1) {

      eta.mu <- X * delta[1]

      X <- cbind(X,rep(0,length(X)))

    }

    else {

      beta <- delta[1:k.beta]

      eta.mu <- X \%*\% beta

    }

    if(k.alpha == 1) {

      eta.phi <- W * delta[k.beta + 1]

      W <- cbind(W,rep(0,length(W)))

    }

    else {

      alpha <- delta[(k.beta + 1) : (k.beta + k.alpha)]

      eta.phi <- W \%*\% alpha

    }

    if(k.gamma == 1) {

      eta.omega <- Z * delta[k.beta + k.alpha + 1]

      Z <- cbind(Z,rep(0,length(Z)))

    }

    else {

      gamma <- delta[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)]

      eta.omega <- Z \%*\% gamma

    }



    mu.i <- t.i*exp(eta.mu)

    k.i <- exp(eta.omega)

    b.i <- exp(eta.phi)

    phi.i <- 1 + b.i

    P0 <- exp(-1/phi.i*mu.i)



# Derivative for beta

for (j in 1:k.beta)

  {

# terms of gradient with 1{Y(i) = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* X[,j]*(P0*(-1/phi.i)*mu.i/(k.i+P0)))

# terms of gradient with 1{Y(i) > 0}

    s2 <- sum(ifelse(Y>0,1,0)* X[,j]*(1+(Y-1)*mu.i/(mu.i+b.i*Y)-1/phi.i*mu.i))

  grad[j] <- s1 + s2

  }





# Derivative for alpha

for (j in 1:k.alpha)

  {

# terms of gradient with 1{Y(i) = 0}

    s1 <- sum(ifelse(Y == 0,1,0)*W[,j]*b.i*(1/phi.i^2*P0*mu.i/(k.i+P0)))

# terms of gradient with 1{Y(i) > 0}

    s2 <- sum(ifelse(Y>0,1,0)*W[,j]*b.i*(Y*(Y-1)/(mu.i + b.i*Y)-Y/phi.i+(mu.i-Y)/phi.i^2))

  grad[k.beta+j] <- s1 + s2

  }





# Derivative for gamma

for (j in 1:k.gamma)

  {

# terms of gradient with 1{Y(i) = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* Z[,j]* k.i * (1/(k.i + P0)))

# terms of gradient independent of Y(i)

    s2 <- -sum(Z[,j]*k.i/(1+k.i))

  grad[k.beta+k.alpha+j] <- s1 + s2

  }



  return(-grad)

  }
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
