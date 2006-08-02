"optimized.run" <-
function(Y,X,W,Z)

{



# improved initial values

# get initial beta, phi & omega

  out <- mle.zigp.full.like(Y, X, Offset = t.i, summary=FALSE)

  beta.start <- out$Coefficients

  phi.first <- out$Dispersion.Parameter

  omega.first <- out$ZI.Parameter



# LM for alpha

  rechte.seite <- rep(log(phi.first - 1),n)

  out<-lm(rechte.seite ~ W-1)

  alpha.start <- out$coefficients



# LM for gamma

  rechte.seite <- rep(log(omega.first)-log(1-omega.first),n)

  out<-lm(rechte.seite ~ Z-1)

  gamma.start <- out$coefficients

  rm(rechte.seite,out)



  start.delta <- c(beta.start,alpha.start,gamma.start)

  rm(beta.start,alpha.start,gamma.start)



  return(start.delta)

}

