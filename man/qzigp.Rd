\name{qzigp}
\alias{qzigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Quantile function of ZIGP distribution }
\description{
'qzigp' calculates the quantiles of the ZIGP distribution.
}
\usage{
qzigp(p, mu, phi, omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ vector of probabilities }
  \item{mu}{ mean }
  \item{phi}{ dispersion parameter }
  \item{omega}{ zero inflation parameter }
}
\value{
Calculates a vector of the same length as of p containing quantiles of the
ZIGP distribution.}
\examples{
p <- seq(0, 1, 0.1)
qzigp(p, 2, 1.5, 0.2)
#[1]  0  0  0  0  0  1  1  2  3  4 75
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
