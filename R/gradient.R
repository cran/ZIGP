"gradient" <-
function(delta)

{

  n <- get("n", pos=globalenv())
  Ysave <- get("Ysave", pos=globalenv())
  k.beta <- get("k.beta", pos=globalenv())
  k.alpha <- get("k.alpha", pos=globalenv())
  k.gamma <- get("k.gamma", pos=globalenv())
  t.i <- get("t.i", pos=globalenv())
  
  grad <- double(k.beta + k.alpha + k.gamma)

  eta.mu <- double(n)

  eta.phi <- double(n)

  eta.omega <- double(n)

  s1 <- double(1)

  s2 <- double(1)



    if(k.beta == 1) {

      eta.mu <- Xsave * delta[1]

    }

    else {

      beta <- delta[1:k.beta]

      eta.mu <- Xsave %*% beta

    }

  if(is.null(Wsave)==FALSE){

    if(k.alpha == 1) {

      eta.phi <- Wsave * delta[k.beta + 1]

    }

    else {

      alpha <- delta[(k.beta + 1) : (k.beta + k.alpha)]

      eta.phi <- Wsave %*% alpha

    }

  }

  if(is.null(Zsave)==FALSE){

    if(k.gamma == 1) {

      eta.omega <- Zsave * delta[k.beta + k.alpha + 1]

    }

    else {

      gamma <- delta[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)]

      eta.omega <- Zsave %*% gamma

    }

  }



    mu.i <- t.i*exp(eta.mu)

    if(is.null(Zsave)==FALSE) {

      k.i <- exp(eta.omega)

    }
    
    else { k.i <- rep(0,n) }

    if(is.null(Wsave)==FALSE) {

      b.i <- exp(eta.phi)

      phi.i <- 1+b.i

    }

    else { b.i <- rep(0,n)

           phi.i <- rep(1,n) }

    P0 <- exp(-1/phi.i*mu.i)

    P0[P0<10e-100] <- 10e-100

if(k.beta==1){ Xsave <- cbind(Xsave,rep(0,length(Xsave)))  }

if(k.alpha==1){ Wsave <- cbind(Wsave,rep(0,length(Wsave)))  }

if(k.gamma==1){ Zsave <- cbind(Zsave,rep(0,length(Zsave)))  }


# Derivative for beta

for (j in 1:k.beta)

  {

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Ysave == 0,1,0)* Xsave[,j]*(P0*(-1/phi.i)*mu.i/(k.i+P0)))

# terms of gradient with 1_{Y_i > 0}

    s2 <- sum(ifelse(Ysave>0,1,0)* Xsave[,j]*(1+(Ysave-1)*mu.i/(mu.i+b.i*Ysave)-1/phi.i*mu.i))

  grad[j] <- s1 + s2

  }





# Derivative for alpha

if(is.null(Wsave)==FALSE) {

  for (j in 1:k.alpha)

    {

  # terms of gradient with 1_{Y_i = 0}

      s1 <- sum(ifelse(Ysave == 0,1,0)*Wsave[,j]*b.i*(1/phi.i^2*P0*mu.i/(k.i+P0)))

  # terms of gradient with 1_{Y_i > 0}

      s2 <- sum(ifelse(Ysave>0,1,0)*Wsave[,j]*b.i*(Ysave*(Ysave-1)/(mu.i + b.i*Ysave)-Ysave/phi.i+(mu.i-Ysave)/phi.i^2))

    grad[k.beta+j] <- s1 + s2

    }

}





# Derivative for gamma

if(is.null(Zsave)==FALSE) {

  for (j in 1:k.gamma)

    {

  # terms of gradient with 1_{Y_i = 0}

      s1 <- sum(ifelse(Ysave == 0,1,0)* Zsave[,j]* k.i * (1/(k.i + P0)))

  # terms of gradient independent of Y_i

      s2 <- -sum(Zsave[,j]*k.i/(1+k.i))

    grad[k.beta+k.alpha+j] <- s1 + s2

    }

}



  return(-grad)

}

