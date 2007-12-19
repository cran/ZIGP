\name{rzigp}
\alias{rzigp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ZIGP random variable generator }
\description{
'rzigp' generates ZIGP random variables using the inversion method.
}
\usage{
rzigp(n, mu, phi, omega)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ length of output random vector. }
  \item{mu}{ mean }
  \item{phi}{ dispersion parameter }
  \item{omega}{ zero inflation parameter }
}
\value{
Generates a ZIGP random vector of length n.}
\examples{
# set seed for random variable generator
set.seed(2)
rzigp(3, 2, 1.5, 0.2)
#[1] 3 1 0
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
