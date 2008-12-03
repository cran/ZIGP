\name{est.zigp}
\alias{est.zigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting ZIGP(mu(i), phi(i), omega(i)) - Regression Models }
\description{
'est.zigp' is used to fit ZIGP(mu(i), phi(i), omega(i)) - Regression Models.
}
\usage{
est.zigp(Yin, fm.X, fm.W=NULL, fm.Z=NULL, 
          Offset = rep(1, length(Yin)), init = T, tex=F)
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
  \item{tex}{ Should the output be a TeX table? Defaults to regular output. }
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
\seealso{ Explorary data analysis tools eda.od() and eda.zi(). }
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

est.zigp(Yin=Y, fm.X=fm.X, fm.W=fm.W, fm.Z=fm.Z, Offset = insurance.year, 
          init = FALSE)     



# approximate equivalence of Poisson-glm and ZIGP-package results
# glm uses IWLS, ZIGP uses numerical maximization of the log-likelihood
# (time series character of the data is neglected)

data(Seatbelts)
DriversKilled <- as.vector(Seatbelts[,1])            # will be response
kms <- as.vector(Seatbelts[,5]/mean(Seatbelts[,5]))  # will be exposure
PetrolPrice <- as.vector(Seatbelts[,6])              # will be covariate 1
law <- as.vector(Seatbelts[,8])                      # will be covariate 2

fm.X <- DriversKilled ~ PetrolPrice + law
out.glm <- glm(fm.X, family=poisson, offset=log(kms))
summary(out.glm)

fm.X <- ~ PetrolPrice + law
est.zigp(DriversKilled, fm.X = fm.X, NULL, NULL, Offset = kms)

# GP with constant overdispersion
fm.X <- ~ PetrolPrice + law
fm.W <- ~ 1
est.zigp(DriversKilled, fm.X, fm.W, NULL, Offset = kms)

# ZIP with constant zero-inflation
fm.X <- ~ PetrolPrice + law
fm.Z <- ~ 1
est.zigp(DriversKilled, fm.X, NULL, fm.Z, Offset = kms)

# ZIGP with constant overdispersion and constant zero-inflation
fm.X <- ~ PetrolPrice + law
fm.W <- ~ 1
fm.Z <- ~ 1
est.zigp(DriversKilled, fm.X, fm.W, fm.Z, Offset = kms)
# no significant zero-inflation according to the Wald test
# (not surprising since not a single zero outcome in data)

# generate TeX output
est.zigp(DriversKilled, fm.X, fm.W, fm.Z, Offset = kms, tex=TRUE)
}
\references{
Czado, C., Erhardt, V., Min, A., Wagner, S. (2007) 
Zero-inflated generalized Poisson models with regression effects on the mean, 
dispersion and zero-inflation level applied to patent outsourcing rates.
Statistical Modelling 7 (2), 125-153.
}
\keyword{ models }% at least one, from doc/KEYWORDS

