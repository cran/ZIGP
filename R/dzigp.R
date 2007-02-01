dzigp <-
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
