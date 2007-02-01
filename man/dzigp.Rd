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
#[1] 0.201467310 0.144357635 0.093859287 0.058560538 0.035846689 0.021749120
#[7] 0.013146607 0.007938715 0.004796227 0.002901477
 
## The function is currently defined as
function(x, mu = stop("no mu arg"), phi = stop("no phi arg"), 
         omega = stop("no omega arg")){
  # check if parameters are valid
  if(omega < 0) {return("omega has to be in [0,1]!")}
  if(omega > 1) {return("omega has to be in [0,1]!")}

   upper <- max(x)
   p <- double(upper+1)

   #P(X=0)
   p[1] <- omega + (1-omega) * exp(-mu/phi)
   if (upper > 0) {   
      rekursive <- FALSE
      for (i in 1:upper) {
        #P(X=x)
        if (rekursive==FALSE) {
          p[i+1] <- (1-omega)*mu*(mu+(phi-1)*i)^(i-1)/exp(lgamma(i+1))*
                    phi^(-i)*exp(-1/phi*(mu+(phi-1)*i))}
        if (p[i+1]==Inf) { 
          rekursive <- TRUE 
          log.p.alt <- log( (1-omega)*mu*(mu+(phi-1)*(i-1))^(i-2)/
                       exp(lgamma(i-1+1))*
                       phi^(-(i-1))*exp(-1/phi*(mu+(phi-1)*(i-1))))
          }
        if (rekursive==TRUE) {
          log.p <- log( (mu+(i-1)*(phi-1))/(phi*i)*
                   (1+(phi-1)/(mu+(i-1)*(phi-1)))^(i-1)*
                   exp(1/phi-1) ) + log.p.alt
          log.p.alt <- log.p
          p[i+1] <- exp(log.p)
        }
      }
   }
   p2 <- double(length(x))
   for (i in 1:length(x)) {
     p2[i] <- p[x[i]+1]
   }
   return(p2)
}
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
