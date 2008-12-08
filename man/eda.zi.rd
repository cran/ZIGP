\name{eda.zi}
\alias{eda.zi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exploratory data analysis tool for zero-inflation level }
\description{
'eda.zi' performs an exploratory data analysis on the influence of a 
covariate on the zero-inflation design (where the logit-link is assumed).
Thereby, a discretization using scoring classes will be applied and empirical
logits be calculated for each scoring class (see Czado et. al (2007)).
Here, a shift of 1/2 is used to obtain well defined empirical logits even for
0. The dashed line is the empirical logit of 1/(number of scoring classes).
Empirical logits further away from this line indicate high influence on
zero-inflation.
}
\usage{
eda.zi(x, y, numberclasses=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Covariate considered }
  \item{y}{ Response considered }
  \item{numberclasses}{ Number of classes for discretization. Defaults to 5. }
}
\details{
As covariate x, discrete or continuous variables may be considered. Categorical
covariates with two levels are allowed as well.

Notwithstanding the description in Czado et. al (2007), the empirical
logits are now adjusted for individual class sizes.
}
\examples{
data(Seatbelts)
DriversKilled <- as.vector(Seatbelts[,1])           # will be response
kms <- as.vector(Seatbelts[,5]/mean(Seatbelts[,5])) # will be exposure
PetrolPrice <- as.vector(Seatbelts[,6])             # will be covariate 1
law <- as.vector(Seatbelts[,8])                     # will be covariate 2


# artificially create some zeros
DriversKilled[PetrolPrice<0.09] <- 0

eda.zi(x=PetrolPrice, y=DriversKilled)
eda.zi(x=PetrolPrice, y=DriversKilled, numberclasses=200)
eda.zi(x=law, y=DriversKilled)
}
\references{
Czado, C., Erhardt, V., Min, A., Wagner, S. (2007) 
Zero-inflated generalized Poisson models with regression effects on the mean, 
dispersion and zero-inflation level applied to patent outsourcing rates.
Statistical Modelling 7 (2), 125-153.
}
\keyword{ models }% at least one, from doc/KEYWORDS
