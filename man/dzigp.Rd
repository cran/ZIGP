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
x <- 1:10
dzigp(x, 2, 1.5, 0.2)
#[1] 0.201467310 0.144357635 0.093859287 0.058560538 0.03584669 0.02174912
#[7] 0.013146607 0.007938715 0.004796227 0.002901477
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
