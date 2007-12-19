\name{ZIGP-package}
\alias{ZIGP-package}
\alias{ZIGP}
\encoding{latin1}
\docType{package}
\title{
Zero-Inflated Generalized Poisson (ZIGP) Models
}
\description{
Adapts GP, ZIP and ZIGP models by Maximum Likelihood estimation.
Regression is allowed not only on the mean but also on overdispersion
and zero-inflation level. Therefore, three design matrices are allowed.
}
\details{
\tabular{ll}{
Package: \tab ZIGP\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2007-12-19\cr
License: \tab GPL (>= 2)\cr
}
Distribution functions in R notation are 'dzigp', 'pzigp', 'qzigp', 'rzigp'.

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
Czado, C., Erhardt, V., Min, A., Wagner, S. (2007) 
Zero-inflated generalized Poisson models with regression effects on the mean, 
dispersion and zero-inflation level applied to patent outsourcing rates.
Statistical Modelling 7 (2), 125-153.

Masterthesis in German:
Erhardt, Vinzenz. Verallgemeinerte Poisson und Nullenueberschuss- 
Regressionsmodelle mit regressiertem Erwartungswert, Dispersions- und 
Nullenüberschuß-Parameter und eine Anwendung zur Patentmodellierung. 
("http://www-m4.ma.tum.de/Diplarb/"), 2006
}
\keyword{ models }
\keyword{ regression }
