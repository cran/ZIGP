\name{ZIGP-package}
\alias{ZIGP-package}
\alias{ZIGP}
\docType{package}
\title{
Zero-Inflated Generalized Poisson (ZIGP) Models
}
\description{
Adapts ZIGP models by Maximum Likelihood estimation.
Regression is allowed not only on the mean but also on overdispersion
and zero-inflation level. Therefore, three design matrices are used.
}
\details{
\tabular{ll}{
Package: \tab ZIGP\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-08-02\cr
License: \tab GPL\cr
}
The most important function is 'wald.test'. This function can be fed with
design matrices for the mean, overdispersion and zero-inflation level. It 
returns an summary-like overview of the covariates together with significance
statistics referring to the Wald test. It can be used for the sequential 
elimination of non-significant effects.

Other useful functions are 'loglikelihood.zigp', which evaluates the 
loglikelihood function on the given parameter value. Scores can be calculated
using 'gradient', the Fisher Information matrix using 'FM'.
}
\author{
Vinzenz Erhardt

Maintainer: Vinzenz Erhardt <erhardt@ma.tum.de>
}
\references{
Masterthesis in German:
Erhardt, Vinzenz. Verallgemeinerte Poisson und Nullenueberschuss- 
Regressionsmodelle mit regressiertem Erwartungswert, Dispersions- und 
Nullenüberschuß-Parameter und eine Anwendung zur Patentmodellierung. 
("http://www-m4.ma.tum.de/Diplarb/"), 2006
}
\keyword{ models }
\keyword{ regression }
