\name{clarke}
\alias{clarke}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Clarke's test for non-nested model comparison }
\description{
'clarke' suggests the better of two (not necessarily nested) models.
}
\usage{
clarke(model1, model2, alpha=0.05, correction=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model1, model2}{ the output of two model fits obtained by using 'mle.zigp'.}
  \item{alpha}{ significance level, defaults to 0.05.}
  \item{correction}{ boolean, if TRUE (default), the Schwarz correction will be used on the differences of log-likelihoods.}
}
\examples{
data(Seatbelts)
DriversKilled <- as.vector(Seatbelts[,1])           # will be response
kms <- as.vector(Seatbelts[,5]/mean(Seatbelts[,5])) # will be exposure
PetrolPrice <- as.vector(Seatbelts[,6])             # will be covariate 1
law <- as.vector(Seatbelts[,8])                     # will be covariate 2

fm.X.poi <- ~ PetrolPrice + law

fm.X.gp <- ~ PetrolPrice + law
fm.W.gp <- ~ 1

fm.X.zigp <- ~ PetrolPrice + law
fm.W.zigp <- ~ 1
fm.Z.zigp <- ~ 1


poi  <- mle.zigp(Yin=DriversKilled, fm.X=fm.X.poi,  fm.W=NULL,   
                 fm.Z=NULL,  Offset = kms, init = FALSE)
gp   <- mle.zigp(Yin=DriversKilled, fm.X=fm.X.gp,   fm.W=fm.W.gp,   
                 fm.Z=NULL,  Offset = kms, init = FALSE)
zigp <- mle.zigp(Yin=DriversKilled, fm.X=fm.X.zigp, fm.W=fm.W.zigp, 
                 fm.Z=fm.Z.zigp, Offset = kms, init = FALSE)
clarke(poi,gp)
clarke(gp,zigp)
clarke(poi,zigp)
}
\references{
Clarke, Kevin A. (2007). A Simple Distribution-Free Test for Nonnested Model Selection. 
Political Analysis 2007 15(3), 347-363.

Schwarz, G. (1978). Estimating the Dimension of a Model. 
Annals of Statistics 6, 461-464.
}
\keyword{ models }% at least one, from doc/KEYWORDS

