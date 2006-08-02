"rzigp" <-
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

