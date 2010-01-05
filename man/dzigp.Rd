\name{dzigp}
\alias{dzigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ PMF of ZIGP distribution }
\description{
'dzigp' calculates the probability mass function of the ZIGP distribution.
}
\usage{
dzigp(x, mu, phi, omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of discrete points }
  \item{mu}{ mean }
  \item{phi}{ dispersion parameter }
  \item{omega}{ zero inflation parameter }
}
\value{
Calculates a vector of the same length as of x of pmf-values.}
\examples{
x <- 0:10
dzigp(x, 2, 1.5, 0.2)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
