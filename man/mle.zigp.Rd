\name{mle.zigp}
\alias{mle.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimates }
\description{
'mle.zigp' is used to calculate the MLEs of the regression parameters for
mean, overdispersion and zero-inflation.
}
\usage{
mle.zigp(Yin, fm.X, fm.W=NULL, fm.Z=NULL, Offset = rep(1, length(Yin)),
init = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yin}{ response vector of length n. }
  \item{fm.X}{ formula for mean design. }
  \item{fm.W}{ formula for overdispersion design (optional). }
  \item{fm.Z}{ formula for zero inflation design (optional). }
  \item{Offset}{ exposure for individual observation lengths. Defaults to a vector of 1.
        The offset MUST NOT be in 'log' scale. }
  \item{init}{ a logical value indicating whether initial optimization values for
          dispersion are set to -2.5 and values for zero inflation regression 
          parameters are set to -1 (init = F) or are estimated by a 
          ZIGP(mu(i), phi, omega)-model (init = T).  Defaults to 'T'. }
}
\details{
     Constant overdispersion and/or zero-inflation can be modelled using an 
     Intercept design on the corresponding level.
     Setting fm.W to NULL corresponds to modelling a ZIP model.
     Setting fm.Z to NULL corresponds to modelling a GP model.
     Setting fm.W and fm.Z to NULL corresponds to modelling a Poisson GLM.

     For numerical stability it may be very useful to center and standardize
     all non-categorical covariates, i.e. use 'x <- (x-mean(x))/sd(x)'.
}
\examples{
# Number of damages in car insurance.
# (not a good fit, just to illustrate how the software is used)
     
damage <- c(0,1,0,0,0,4,2,0,1,0,1,1,0,2,0,0,1,0,0,1,0,0,0)
insurance.year <- c(1,1.2,0.8,1,2,1,1.1,1,1,1.1,1.2,1.3,0.9,1.4,1,1,1,
1.2,1,1,1,1,1)
drivers.age <- c(25,19,30,48,30,18,19,29,24,54,56,20,38,18,23,58,
47,36,25,28,38,39,42)
# for overdispersion: car brand dummy in {1,2,3}, brand = 1 is reference
brand <- c(1,2,1,3,3,2,2,1,1,3,2,2,1,3,1,3,2,2,1,1,3,3,2)
# abroad: driver has been abroad for longer time (=1)
abroad <- c(0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,1)
Y <- damage
fm.X <- ~ drivers.age
fm.W <- ~ 0 + factor(brand)
fm.Z <- ~ abroad

mle.zigp(Yin=Y, fm.X=fm.X, fm.W=fm.W, fm.Z=fm.Z, Offset = insurance.year, init = FALSE)

}
\references{
Czado, C., Erhardt, V., Min, A., Wagner, S. (2007) 
Zero-inflated generalized Poisson models with regression effects on the mean, 
dispersion and zero-inflation level applied to patent outsourcing rates.
Statistical Modelling 7 (2), 125-153.
}
\keyword{ models }% at least one, from doc/KEYWORDS
