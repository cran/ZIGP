"response.zigp" <-
function(X, W, Z, beta, alpha, gamma)

{

# get number of observations from X, W and Z

  if(is.matrix(X)) {

    nx <- dim(X)[1]

  }

  else {

    nx <- length(X)

  }

  if(is.matrix(W)) {

    nw <- dim(W)[1]

  }

  else {

    nw <- length(W)

  }

  if(is.matrix(Z)) {

    nz <- dim(Z)[1]

  }

  else {

    nz <- length(Z)

  }



  if (nx==nz&nx==nw) {

  n <- nx

  mu <- double(n)

  phi <- double(n)

  omega <- double(n)

  Y <- double(n)



  if(is.matrix(X)) {

    mu <- exp(X%*%beta)

  }

  else {

    mu <- exp(X*beta)

  }

  if(is.matrix(W)) {

    phi <- 1 + exp(W%*%alpha)

  }

  else {

    phi <- 1 + exp(W*alpha)

  }

  if(is.matrix(Z)) {

    omega <- exp(Z%*%gamma)/(1+exp(Z%*%gamma))

  }

  else {

    omega <- exp(Z*gamma)/(1+exp(Z*gamma))

  }

# create ZIGP response for all i

  for(i in 1:n) {

    Y[i] <- rzigp(1, mu[i], phi[i], omega[i])

  }

  return(Y)

  }

  else {

  print("X, W and Z do not have the same numbers of rows!")

  }

}

