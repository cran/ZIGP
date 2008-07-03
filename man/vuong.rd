\name{vuong}
\alias{vuong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Vuong test for non-nested model comparison }
\description{
'vuong' suggests the better of two (not necessarily nested) models according to
Vuong's statistic.
}
\usage{
vuong(model1, model2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model1, model2}{ the output of two model fits obtained by using 'mle.zigp'.
}
}
\examples{
## Number of damages in car insurance.

damage <- c(0,1,0,0,0,4,2,0,1,0,1,1,0,2,0,0,1,0,0,1,0,0,0)
Intercept <- rep(1,length(damage))
insurance.year <- c(1,1.2,0.8,1,2,1,1.1,1,1,1.1,1.2,1.3,0.9,1.4,1,1,1,1.2,
1,1,1,1,1)
drivers.age <- c(25,19,30,48,30,18,19,29,24,54,56,20,38,18,23,58,
47,36,25,28,38,39,42)
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

out1 <- mle.zigp(Yin=Y, Xin=X, Win=W, Zin=Z, Offset = insurance.year, init = FALSE)
out2 <- mle.zigp(Yin=Y, Xin=X, Win=NULL, Zin=NULL, Offset = insurance.year, init = FALSE)
vuong(out1,out2)
#[1] "nu = -0.836"
#[1] "None of the models is favoured."
#[1] "P-value: 0.40"
}
\keyword{ models }% at least one, from doc/KEYWORDS

