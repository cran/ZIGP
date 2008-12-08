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
x <- 1:10
pzigp(x, 2, 1.5, 0.2)
#[1] 0.6123450 0.7567027 0.8505619 0.909123 0.944969 0.966718 0.979865
#[8] 0.9878036 0.9925998 0.9955013
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
