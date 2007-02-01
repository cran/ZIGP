pzigp <-
function(x, mu = stop("no mu arg"), phi = stop("no phi arg"), 
         omega = stop("no omega arg")){
  # check if parameters are valid
  if(omega < 0) {return("omega has to be in [0,1]!")}
  if(omega > 1) {return("omega has to be in [0,1]!")}

   upper <- max(x)
   s <- double(upper+1)

   #P(X=0)
   p <- omega + (1-omega) * exp(-mu/phi)
   s[1] <- p
   if (upper > 0) {   
      rekursive <- FALSE
      for (i in 1:upper) {
        #P(X=x)
        if (rekursive==FALSE) {
          p <- (1-omega)*mu*(mu+(phi-1)*i)^(i-1)/exp(lgamma(i+1))*
               phi^(-i)*exp(-1/phi*(mu+(phi-1)*i))}
        if (p==Inf) { 
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
          p <- exp(log.p)
        }
        s[i+1] <- s[i] + p
      }
   }
   s2 <- double(length(x))
   for (i in 1:length(x)) {
     s2[i] <- s[x[i]+1]
   }
   return(s2)
}
