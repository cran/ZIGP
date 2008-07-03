"gradient" <-
function(delta)

{

  n <- get("n")
  k.beta <- get("k.beta")
  k.alpha <- get("k.alpha")
  k.gamma <- get("k.gamma")
  Y <- get("Y")
  t.i <- get("t.i")
  
  grad <- double(k.beta + k.alpha + k.gamma)

  eta.mu <- double(n)

  eta.phi <- double(n)

  eta.omega <- double(n)

  s1 <- double(1)

  s2 <- double(1)



    if(k.beta == 1) {

      eta.mu <- X * delta[1]

    }

    else {

      beta <- delta[1:k.beta]

      eta.mu <- X %*% beta

    }

  if(is.null(W)==FALSE){

    if(k.alpha == 1) {

      eta.phi <- W * delta[k.beta + 1]

    }

    else {

      alpha <- delta[(k.beta + 1) : (k.beta + k.alpha)]

      eta.phi <- W %*% alpha

    }

  }

  if(is.null(Z)==FALSE){

    if(k.gamma == 1) {

      eta.omega <- Z * delta[k.beta + k.alpha + 1]

    }

    else {

      gamma <- delta[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)]

      eta.omega <- Z %*% gamma

    }

  }



    mu.i <- t.i*exp(eta.mu)

    if(is.null(Z)==FALSE) {

      k.i <- exp(eta.omega)

    }
    
    else { k.i <- rep(0,n) }

    if(is.null(W)==FALSE) {

      b.i <- exp(eta.phi)

      phi.i <- 1+b.i

    }

    else { b.i <- rep(0,n)

           phi.i <- rep(1,n) }

    P0 <- exp(-1/phi.i*mu.i)

if(k.beta==1){ X <- cbind(X,rep(0,length(X)))  }

if(k.alpha==1){ W <- cbind(W,rep(0,length(W)))  }

if(k.gamma==1){ Z <- cbind(Z,rep(0,length(Z)))  }


# Derivative for beta

for (j in 1:k.beta)

  {

# terms of gradient with 1_{Y_i = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* X[,j]*(P0*(-1/phi.i)*mu.i/(k.i+P0)))

# terms of gradient with 1_{Y_i > 0}

    s2 <- sum(ifelse(Y>0,1,0)* X[,j]*(1+(Y-1)*mu.i/(mu.i+b.i*Y)-1/phi.i*mu.i))

  grad[j] <- s1 + s2

  }





# Derivative for alpha

if(is.null(W)==FALSE) {

  for (j in 1:k.alpha)

    {

  # terms of gradient with 1_{Y_i = 0}

      s1 <- sum(ifelse(Y == 0,1,0)*W[,j]*b.i*(1/phi.i^2*P0*mu.i/(k.i+P0)))

  # terms of gradient with 1_{Y_i > 0}

      s2 <- sum(ifelse(Y>0,1,0)*W[,j]*b.i*(Y*(Y-1)/(mu.i + b.i*Y)-Y/phi.i+(mu.i-Y)/phi.i^2))

    grad[k.beta+j] <- s1 + s2

    }

}





# Derivative for gamma

if(is.null(Z)==FALSE) {

  for (j in 1:k.gamma)

    {

  # terms of gradient with 1_{Y_i = 0}

      s1 <- sum(ifelse(Y == 0,1,0)* Z[,j]* k.i * (1/(k.i + P0)))

  # terms of gradient independent of Y_i

      s2 <- -sum(Z[,j]*k.i/(1+k.i))

    grad[k.beta+k.alpha+j] <- s1 + s2

    }

}



  return(-grad)

}

