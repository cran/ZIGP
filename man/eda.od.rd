\name{eda.od}
\alias{eda.od}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exploratory data analysis tool for overdispersion level }
\description{
'eda.od' performs an impact study on the influence of a covariate on the
overdispersion design (where the shifted log-link is assumed). Thereby, a
discretation using scoring classes will be applied and the overdispersion
function be calculated for each scoring class (see Czado et. al (2007)).
}
\usage{
eda.od(x, y, Offset=rep(1,length(y)), numberclasses=5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Covariate considered }
  \item{y}{ Response considered }
  \item{Offset}{ Exposure for individual observation lengths. Defaults to a vector of 1.
          The offset MUST NOT be in 'log' scale. }
  \item{numberclasses}{ Number of classes for discretization. Defaults to 5. }
}
\details{
As covariate x, discrete or continuous variables may be considered. Categorical
covariates only make sense if they have only two levels.
}
\examples{
data(Seatbelts)
DriversKilled <- as.vector(Seatbelts[,1])           # will be response
kms <- as.vector(Seatbelts[,5]/mean(Seatbelts[,5])) # will be exposure
PetrolPrice <- as.vector(Seatbelts[,6])             # will be covariate 1
law <- as.vector(Seatbelts[,8])                     # will be covariate 2

eda.od(x=PetrolPrice, y=DriversKilled, Offset=kms)
eda.od(x=PetrolPrice, y=DriversKilled, Offset=kms, numberclasses=20)
eda.od(x=law, y=DriversKilled, Offset=kms)
}
\references{
Czado, C., Erhardt, V., Min, A., Wagner, S. (2007) 
Zero-inflated generalized Poisson models with regression effects on the mean, 
dispersion and zero-inflation level applied to patent outsourcing rates.
Statistical Modelling 7 (2), 125-153.
}
\keyword{ models }% at least one, from doc/KEYWORDS

