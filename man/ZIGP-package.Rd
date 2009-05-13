\name{ZIGP-package}
\alias{ZIGP-package}
\alias{ZIGP}
\encoding{latin1}
\docType{package}
\title{
Zero-Inflated Generalized Poisson (ZIGP) Models
}
\description{
Fits GP, ZIP and ZIGP models by Maximum Likelihood estimation.
Regression is allowed not only on the mean but also on overdispersion
and zero-inflation level.
}
\details{
\tabular{ll}{
Package: \tab ZIGP\cr
Type: \tab Package\cr
Version: \tab 3.5\cr
Date: \tab 2009-05-13\cr
License: \tab GPL (>= 3)\cr
}
Distribution functions in R notation are 'dzigp', 'pzigp', 'qzigp', 'rzigp'.

The main function is 'est.zigp'. This function can be fed with
design formulas for the mean, overdispersion and zero-inflation level. It 
returns a summary-like overview of the covariates together with significance
statistics referring to the Wald test. It can be used for sequential
elimination of non-significant effects. Function 'mle.zigp' returns estimates
of regression coefficients, AIC etc. for further use.

Other useful functions are 'loglikelihood.zigp', which evaluates the 
loglikelihood function on the given parameter value. Scores can be calculated
using 'gradient', the Fisher Information matrix using 'FM'.

Tools for an exploratory data analysis for the overdispersion level and zero-
inflation level are given in 'eda.od' and 'eda.zi', respectively.

Nonnested model comparison (also for the Negative Binomial distribution) can
be facilitated using a test proposed by Vuong which is implemented in function
'vuong' or by a test proposed by Clarke using function 'clarke'.
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
("http://www-m4.ma.tum.de/Diplarb/"), 2006.

Vuong, Q.H. (1989). Likelihood Ratio tests for model selection and nonnested
hypotheses.
Econometrica 57(2), 307-333.

Clarke, Kevin A. (2007). A Simple Distribution-Free Test for Nonnested Model Selection. 
Political Analysis 2007 15(3), 347-363.

Schwarz, G. (1978). Estimating the Dimension of a Model. 
Annals of Statistics 6, 461-464.
}
\keyword{ models }
\keyword{ regression }
