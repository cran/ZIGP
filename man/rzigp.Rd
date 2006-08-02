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

## The function is currently defined as
function(n, mu = stop("no mu arg"), phi = stop("no phi arg"), omega = stop("no omega arg"))

{

# check if parameters are valid

if(omega < 0) {return("omega has to be in [0,1]!")}

if(omega > 1) {return("omega has to be in [0,1]!")}



# inversion method

    x <- double(n)

    for(i in 1:n) {

      y <- integer(1)

      c <- double(1)

      #P(y=0)

      p <- omega + (1-omega) * exp( - mu/phi)

      s <- p

      u <- runif(1, 0, 1)

      if (u > s){

        #P(y=1)

        p <- (1-omega)*mu/phi*exp(-1/phi*(mu+phi-1))

        s <- s + p

        y <- 1

        while (u > s) {

          y <- y + 1

          p <- (mu + (phi-1)*y)^(y-1) / (mu + (phi-1)*(y-1))^(y-2) / (y * phi) *

          exp((1-phi)/phi) * p

          s <- s + p

        }

      }

      x[i] <- y

    }

    return(x)

  }
}
\keyword{ datagen }% at least one, from doc/KEYWORDS
