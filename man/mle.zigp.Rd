\name{mle.zigp}
\alias{mle.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimates }
\description{
'mle.zigp' is used to calculate the MLEs of the regression parameters for
mean, overdispersion and zero-inflation.
}
\usage{
mle.zigp(Yin, Xin, Win=NULL, Zin=NULL, Offset = rep(1, length(Yin)),
init = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yin}{ response vector of length n. }
  \item{Xin}{ design matrix of dim (n x p) for mean modelling. }
  \item{Win}{ design matrix of dim (n x r) for overdispersion modelling. }
  \item{Zin}{ design matrix of dim (n x q) for zero inflation modelling. }
  \item{Offset}{ exposure for individual observation lengths. Defaults to a vector of 1.
        The offset MUST NOT be in 'log' scale. }
  \item{init}{ a logical value indicating whether initial optimization values for
          dispersion are set to -2.5 and values for zero inflation regression 
          parameters are set to -1 (init = F) or are estimated by a 
          ZIGP(mu(i), phi, omega)-model (init = T).  Defaults to 'T'. }
}
\examples{
## Number of damages in car insurance.

damage <- c(0,1,0,0,0,4,2,0,1,0,1,1,0,2,0,0,1,0,0,1,0,0,0)
Intercept <- rep(1,length(damage))
insurance.year <- c(1,1.2,0.8,1,2,1,1.1,1,1,1.1,1.2,1.3,0.9,1.4,1,1,1,1.2,
1,1,1,1,1)
drivers.age <- c(25,19,30,48,30,18,19,29,24,54,56,20,38,18,23,58,
47,36,25,28,38,39,42)
n <- length(damage)
# for overdispersion: car brand dummy in {1,2,3}, brand = 1 is reference
brand <- c(1,2,1,3,3,2,2,1,1,3,2,2,1,3,1,3,2,2,1,1,3,3,2)
brand2 <- ifelse(brand==2,1,0)
brand3 <- ifelse(brand==3,1,0)
W <- cbind(brand2,brand3)
# abroad: driver has been abroad for longer time (=1)
abroad <- c(0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,1)
Y <- damage
X <- cbind(Intercept, drivers.age)
Z <- cbind(abroad)

mle.zigp(Yin=Y, Xin=X, Win=W, Zin=Z, Offset = insurance.year, init = FALSE)

#[1]   Range for ZI-Parameters:   0.2491062 0.5
#[2]   Range of Dispersion Pars:  1.000176 2.189325
#[3]   Coefficients for mu:       1.471478 -0.05075418
#[4]   Coefficients for phi:      -8.646371 0.1733860
#[5]   Coefficients for omega:    -1.103385
#[6]   Pearson Chi Squared:       15.09779
#[7]   Range of mu:               0.2294054 2.445806
#[8]   Message:                   "NULL"
#[9]   AIC:                       56.88305
}
\keyword{ models }% at least one, from doc/KEYWORDS
