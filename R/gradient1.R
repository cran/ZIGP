"gradient1" <-
function(delta)

{

  grad <- double(2+k)

  eta <- double(n)

  s1 <- double(1)

  s2 <- double(1)



    if(k == 1) {

      eta <- X * delta[3]

      X <- cbind(X,rep(0,length(X)))

    }

    else {

      beta <- delta[3:(k+2)]

      eta <- X %*% beta

    }



    mu.i <- t.i*exp(eta)

    phi <- 1+exp(delta[1])

    omega <- exp(delta[2])/(1+exp(delta[2]))

    P0 <- exp(-1/phi*mu.i)



# Derivative for phi

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* (1/phi^2*P0*mu.i/(omega/(1-omega)+P0)))

# terms of gradient with 1_{Y_i > 0}

    s2 <- sum(ifelse(Y>0,1,0)* (Y*(Y-1)/(mu.i+(phi-1)*Y)-Y/phi+(mu.i-Y)/phi^2))

  grad[1] <- s1 + s2



# Derivative for omega

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* 1/(((1-omega)^2)*(omega/(1-omega) + P0)))

# terms of gradient independent of Y_i

    s2 <- -n/(1-omega)

  grad[2] <- s1 + s2



# Derivative for beta

for (j in 1:k)

  {

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* X[,j]*(P0*(-1/phi)*mu.i/(omega/(1-omega)+P0)))

# terms of gradient with 1_{Y_i > 0}

    s2 <- sum(ifelse(Y>0,1,0)* X[,j]*(1+(Y-1)*mu.i/(mu.i+(phi-1)*Y)-1/phi*mu.i))

  grad[2+j] <- s1 + s2

  }



  return(-grad)

}

