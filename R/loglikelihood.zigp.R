"loglikelihood.zigp" <-
function(delta)

{

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

    

    if(k.alpha == 1) {

      eta.phi <- W * delta[k.beta + 1]

    }

    else {

      alpha <- delta[(k.beta + 1) : (k.beta + k.alpha)]

      eta.phi <- W %*% alpha

    }

    

    if(k.gamma == 1) {

      eta.omega <- Z * delta[k.beta + k.alpha + 1]

    }

    else {

      gamma <- delta[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)]

      eta.omega <- Z %*% gamma

    }





    mu.i <- t.i*exp(eta.mu)

    b.i <- exp(eta.phi)

    phi.i <- 1 + b.i

    k.i <- exp(eta.omega)

    temp1 <- (mu.i + b.i*Y)



# terms of log likelihood with 1_{Y_i = 0}

    s1 <- sum(ifelse(Y == 0,1,0)* (log(k.i+exp(-1/phi.i*mu.i)) -

          log(1+k.i)))

# terms of log likelihood with 1_{Y_i > 0} (excluding factorial)

    s2 <- sum(ifelse(Y>0,1,0)*(-log(1+k.i)+log(t.i)+eta.mu+(Y-1)*log(temp1)-

          Y*log(phi.i)-1/phi.i*temp1) - lgamma(Y+1))



  l <- s1 + s2

  return( - l)

}
