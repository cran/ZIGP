"FM" <-
function(beta, alpha, gamma, X, W, Z, Offset = NULL)

{



k.beta <- length(beta)

k.alpha <- length(alpha)

k.gamma <- length(gamma)



FM <- array(0,c(k.beta+k.alpha+k.gamma,k.beta+k.alpha+k.gamma))



if(is.matrix(X)){ eta.mu <- X%*%beta }

else{eta.mu <- X*beta}

if(is.null(W)==FALSE) {

  if(k.alpha == 1) {

    eta.phi <- W * alpha

  }

  else {

    eta.phi <- W %*% alpha

  }

}

if(is.null(Z)==FALSE) {

  if(k.gamma == 1) {

    eta.omega <- Z * gamma

  }

  else {

    eta.omega <- Z %*% gamma

  }

}



if (is.null(Offset)) {mu  <- exp(eta.mu)}

else {t.i <<- Offset

     mu  <- t.i*exp(eta.mu)}

if(is.null(W)==FALSE) {

  b <- exp(eta.phi)

  phi <- 1+b

}

else { b <- rep(0,n)

       phi <- rep(1,n) }
       
if(is.null(Z)==FALSE) {

  k <- exp(eta.omega)
  
}

else { k <- rep(0,n) }

P0 <- exp(-1/phi*mu)

# a_i = E[I_{Y_i = 0}]

a <- (k+P0)/(1+k)



temp <- double(1)



if(k.beta==1){ X <- cbind(X,rep(0,length(X)))  }

if(k.alpha==1){ W <- cbind(W,rep(0,length(W)))  }

if(k.gamma==1){ Z <- cbind(Z,rep(0,length(Z)))  }



# the second derivative w.r.t. "beta"

for (i in 1:k.beta)

{

  for (j in 1:k.beta)

  {

  if(is.null(Z)==FALSE) {

    temp <- sum( X[,i]*X[,j]*mu* ( a*(-1/phi*P0^2+(mu-phi)/(phi^2)*P0*k)/

            ((k+P0)^2) + (b*mu)/(phi^2*(mu-2+2*phi)*(1+k)) - (1-a)/phi ) )
            
  }
  
  else {
  
    temp <- sum( X[,i]*X[,j]*mu* ( a*(-1/phi) +

            (b*mu)/(phi^2*(mu-2+2*phi)*(1+k)) - (1-a)/phi ) )

  }

  FM[i,j] <- -temp

  }

}



# the second derivative w.r.t. "alpha"

if(is.null(W)==FALSE) {

for (i in 1:k.alpha)

{

  for (j in 1:k.alpha)

  {

  temp <- sum( W[,i]*W[,j]*b*( a*P0*mu/((k+P0)^2)*(mu*b*k/phi^4+(1/phi^2-2*b/phi^3)*(k+P0))+

          mu^2/(phi^2*(mu-2+2*phi)*(1+k)) + 2*mu/(1+k)*(-1/phi^2+b/phi^3) +

          (1-a)/(phi^3)*(mu*phi-2*mu*b) ) )

  FM[k.beta+i,k.beta+j] <- -temp

  }

}

}



# the second derivative w.r.t. "gamma"

if(is.null(Z)==FALSE) {

for (i in 1:k.gamma)

{

  for (j in 1:k.gamma)

  {

  temp <- sum( Z[,i]*Z[,j]*k* (a*P0/((k+P0)^2)-1/((1+k)^2) ) )

  FM[k.beta+k.alpha+i,k.beta+k.alpha+j] <- -temp

  }

}

}



# the mixed derivative w.r.t. "beta" and "alpha"

if(is.null(W)==FALSE) {

for (i in 1:k.beta)

{

  for (j in 1:k.alpha)

  {

temp <- sum( X[,i]*W[,j]*mu*b* ( a*(-1/phi^3*P0*mu*k+1/(phi^2)*P0*(k+P0))/((k+P0)^2)-

        mu/(phi^2*(mu-2+2*phi)*(1+k))+(1-a)/(phi^2) ) )

FM[k.beta+j,i] <- -temp

FM[i,k.beta+j] <- -temp

  }

}

}



# the mixed derivative w.r.t. "beta" and "gamma"

if(is.null(Z)==FALSE) {

for (i in 1:k.beta)

{

  for (j in 1:k.gamma)

  {

  temp <- sum( X[,i]*Z[,j]*a*P0*mu*k/(phi*(k+P0)^2) )

  FM[k.beta+k.alpha+j,i] <- -temp

  FM[i,k.beta+k.alpha+j] <- -temp

  }

}

}



# the mixed derivative w.r.t. "alpha" and "gamma"

if(is.null(W)==FALSE&is.null(Z)==FALSE) {

for (i in 1:k.gamma)

{

  for (j in 1:k.alpha)

  {

  temp <- sum( W[,j]*Z[,i]*b*a* (-1/phi^2 * P0 * mu * k)/((k+P0)^2) )

  FM[k.beta+k.alpha+i,k.beta+j] <- -temp

  FM[k.beta+j,k.beta+k.alpha+i] <- -temp

  }

}

}



return(FM)

}

