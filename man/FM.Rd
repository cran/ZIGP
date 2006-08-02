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
beta  <- c(5,-2)
alpha <- 3.4
gamma <- -10
FM(beta, alpha, gamma, X, W, Z)
#              [,1]          [,2]         [,3]          [,4]
#[1,]  1.469180e+03  2.412237e+03 6.249334e-03 -8.400641e-11
#[2,]  2.412237e+03  3.965799e+03 1.040260e-02 -3.360257e-10
#[3,]  6.249334e-03  1.040260e-02 2.101615e-03  2.520099e-10
#[4,] -8.400641e-11 -3.360257e-10 2.520099e-10  9.079162e-05
     
## The function is currently defined as
function(beta, alpha, gamma, X, W, Z, Offset = NULL)

{



k.beta <- length(beta)

k.alpha <- length(alpha)

k.gamma <- length(gamma)



FM <- array(0,c(k.beta+k.alpha+k.gamma,k.beta+k.alpha+k.gamma))



if(is.matrix(X)){ eta.mu <- X\%*\%beta }

else{eta.mu <- X*beta}

if(is.matrix(W)){ eta.phi <- W\%*\%alpha }

else{eta.phi <- W*alpha}

if(is.matrix(Z)){ eta.gamma <- Z\%*\%gamma }

else{eta.gamma <- Z*gamma}



if (is.null(Offset)) {mu  <- exp(eta.mu)}

else {t.i <<- Offset

     mu  <- t.i*exp(eta.mu)}

b <- exp(eta.phi)

phi <- 1+b

k <- exp(eta.gamma)

P0 <- exp(-1/phi*mu)

# a(i) = E[I{Y(i) = 0}]

a <- (k+P0)/(1+k)



temp <- double(1)



if(k.beta==1){ X <- cbind(X,rep(0,length(X)))  }

if(k.alpha==1){ W <- cbind(W,rep(0,length(W)))  }

if(k.gamma==1){ Z <- cbind(Z,rep(0,length(Z)))  }



# the second derivative w.r.t. "beta"

for (i in 1:k.beta)

{

  for (j in 1:k.beta)

  {

  temp <- sum( X[,i]*X[,j]*mu* ( a*(-1/phi*P0^2+(mu-phi)/(phi^2)*P0*k)/

          ((k+P0)^2) + (b*mu)/(phi^2*(mu-2+2*phi)*(1+k)) - (1-a)/phi ) )

  FM[i,j] <- -temp

  }

  }



# the second derivative w.r.t. "alpha"

for (i in 1:k.alpha)

{

  for (j in 1:k.alpha)

  {

  temp <- sum( W[,i]*W[,j]*b*( a*P0*mu/((k+P0)^2)*(mu*b*k/phi^4+(1/phi^2-2*b/phi^3)*(k+P0))+

          mu^2/(phi^2*(mu-2+2*phi)*(1+k)) + 2*mu/(1+k)*(-1/phi^2+b/phi^3) +

          (1-a)/(phi^3)*(mu*phi-2*mu*b) ) )

  FM[k.beta+i,k.beta+j] <- -temp

  }

  }



# the second derivative w.r.t. "gamma"

for (i in 1:k.gamma)

{

  for (j in 1:k.gamma)

  {

  temp <- sum( Z[,i]*Z[,j]*k* (a*P0/((k+P0)^2)-1/((1+k)^2) ) )

  FM[k.beta+k.alpha+i,k.beta+k.alpha+j] <- -temp

  }

  }



# the mixed derivative w.r.t. "beta" and "alpha"

for (i in 1:k.beta)

{

  for (j in 1:k.alpha)

  {

temp <- sum( X[,i]*W[,j]*mu*b* ( a*(-1/phi^3*P0*mu*k+1/(phi^2)*P0*(k+P0))/((k+P0)^2)-

        mu/(phi^2*(mu-2+2*phi)*(1+k))+(1-a)/(phi^2) ) )

FM[k.beta+j,i] <- -temp

FM[i,k.beta+j] <- -temp

  }

  }



# the mixed derivative w.r.t. "beta" and "gamma"

for (i in 1:k.beta)

{

  for (j in 1:k.gamma)

  {

  temp <- sum( X[,i]*Z[,j]*a*P0*mu*k/(phi*(k+P0)^2) )

  FM[k.beta+k.alpha+j,i] <- -temp

  FM[i,k.beta+k.alpha+j] <- -temp

  }

  }



# the mixed derivative w.r.t. "alpha" and "gamma"

for (i in 1:k.gamma)

{

  for (j in 1:k.alpha)

  {

  temp <- sum( W[,j]*Z[,i]*b*a* (-1/phi^2 * P0 * mu * k)/((k+P0)^2) )

  FM[k.beta+k.alpha+i,k.beta+j] <- -temp

  FM[k.beta+j,k.beta+k.alpha+i] <- -temp

  }

  }



return(FM)

  }
}
\keyword{ distribution }% at least one, from doc/KEYWORDS

