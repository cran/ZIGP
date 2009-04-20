"gradient1" <-
function(delta)

{

  k <- get("ksave", pos=globalenv())
  n <- get("n", pos=globalenv())
  t.i <- get("t.i", pos=globalenv())
  Ysave <- get("Ysave", pos=globalenv())
  grad <- double(2+k)

  eta <- double(n)

  s1 <- double(1)

  s2 <- double(1)



    if(k == 1) {

      eta <- Xsave * delta[3]

      Xsave <- cbind(Xsave,rep(0,length(Xsave)))

    }

    else {

      beta <- delta[3:(k+2)]

      eta <- Xsave %*% beta

    }



    mu.i <- t.i*exp(eta)

    phi <- 1+exp(delta[1])

    omega <- exp(delta[2])/(1+exp(delta[2]))

    P0 <- exp(-1/phi*mu.i)

    P0[P0<10e-100] <- 10e-100



# Derivative for phi

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Ysave == 0,1,0)* (1/phi^2*P0*mu.i/(omega/(1-omega)+P0)))

# terms of gradient with 1_{Y_i > 0}

    s2 <- sum(ifelse(Ysave>0,1,0)* (Ysave*(Ysave-1)/(mu.i+(phi-1)*Ysave)-Ysave/phi+(mu.i-Ysave)/phi^2))

  grad[1] <- s1 + s2



# Derivative for omega

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Ysave == 0,1,0)* 1/(((1-omega)^2)*(omega/(1-omega) + P0)))

# terms of gradient independent of Y_i

    s2 <- -n/(1-omega)

  grad[2] <- s1 + s2



# Derivative for beta

for (j in 1:k)

  {

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Ysave == 0,1,0)* Xsave[,j]*(P0*(-1/phi)*mu.i/(omega/(1-omega)+P0)))

# terms of gradient with 1_{Y_i > 0}

    s2 <- sum(ifelse(Ysave>0,1,0)* Xsave[,j]*(1+(Ysave-1)*mu.i/(mu.i+(phi-1)*Ysave)-1/phi*mu.i))

  grad[2+j] <- s1 + s2

  }



  return(-grad)

}

