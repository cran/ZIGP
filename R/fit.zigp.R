"fit.zigp" <-
function(delta)

{

    n <- get("n", pos=globalenv())
    k.beta <- get("k.beta", pos=globalenv())
    k.alpha <- get("k.alpha", pos=globalenv())
    k.gamma <- get("k.gamma", pos=globalenv())
    X <- get("X", pos=globalenv())
    W <- get("W", pos=globalenv())
    Z <- get("Z", pos=globalenv())
    t.i <- get("t.i", pos=globalenv())
    
    eta.mu <- double(n)

    eta.phi <- double(n)

    eta.omega <- double(n)

    mu <- double(n)

    phi <- double(n)

    omega <- double(n)

    fit <- double(n)



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



# computation of mu, phi, omega and mean

    mu <- t.i*exp(eta.mu)

    if(is.null(W)==FALSE){ phi <- 1 + exp(eta.phi) }
    
    else { phi <- 1 }

    if(is.null(W)==FALSE){ omega <- exp(eta.omega)/(1+exp(eta.omega)) }
    
    else { omega <- 0 }

    fit <- (1 - omega) * mu



        return(list(fit = fit, mu = mu, phi = phi, omega = omega))

}

