\name{pzigp}
\alias{pzigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Distribution function of ZIGP distribution }
\description{
'pzigp' calculates the distribution function of the ZIGP distribution.
}
\usage{
pzigp(x, mu, phi, omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of discrete points }
  \item{mu}{ mean }
  \item{phi}{ dispersion parameter }
  \item{omega}{ zero inflation parameter }
}
\value{
Calculates a vector of the same length as of x evaluating the ZIGP
distribution function at x.}
\examples{
x <- 0:10
pzigp(x, 2, 1.5, 0.2)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
