"mle.zigp" <-
function(Yin, Xin, Win, Zin, Offset = rep(1,length(Yin)), summary = TRUE, init = FALSE)

{

  Y <<- Yin

  X <<- Xin

  W <<- Win

  Z <<- Zin

  if(is.matrix(X)) {

    n <- dim(X)[1]

    k.beta <<- dim(X)[2]

  }

  else {

    n <- length(X)

    k.beta <<- 1

  }

  if(is.matrix(W)) {

    k.alpha <<- dim(W)[2]

    nw <<- dim(W)[1]

  }

  else {

    k.alpha <<- 1

    nw <<- length(W)

  }

  if(is.matrix(Z)) {

    k.gamma <<- dim(Z)[2]

    nz <<- dim(Z)[1]

  }

  else {

    k.gamma <<- 1

    nz <<- length(Z)

  }



  if(nz==n & nw==n){

  t.i <<- Offset



  if(init){   #optimized initial values

    start.delta <- optimized.run(Y,X,W,Z)



    # maximization of log likelihood

    start <- proc.time()[2]

    opt <- optim(par = start.delta, fn = loglikelihood.zigp, gr=gradient,

    method = "BFGS")

    #  method = "CG")

  }

  else

  {

    # simple initial values

    alpha <- rep(-2.5,k.alpha)

    gamma <- rep(-2.5,k.gamma)

    beta <- summary(glm(Y ~ offset(log(t.i)) + 1 + X, family = poisson(link=log)

                    ))$coefficients[, 1]

    start.delta <- c(beta,alpha,gamma)

    # maximization of log likelihood w.r.t. delta s.t. above limits hold

    start <- proc.time()[2]

    opt <- optim(par = start.delta, fn = loglikelihood.zigp, gr=gradient,

    method = "BFGS")

    #  method = "CG")

  }



  delta <- opt$par

  it <- opt$counts

  loglikelihood <-  - opt$value

  if(is.null(opt$message)) {message <- "NULL"}

  else { message <- opt$message }

  AIC <- -2 * loglikelihood + 2 * (k.beta + k.alpha + k.gamma)

  time.needed <- proc.time()[2] - start



# results

  coef <- double(k.beta + k.alpha + k.gamma)

  coef <- delta

  beta <- coef[1 : k.beta]

  alpha <- coef[(k.beta + 1) : (k.beta + k.alpha)]

  gamma <- coef[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)]



  fit <- fit.zigp(delta)

  res <- Y - fit$fit

  chsq <- sum((Y - fit$fit) / sqrt((1-fit$omega)*fit$phi*fit$fit+fit$omega*(1-fit$omega)*fit$fit^2))



  range.mu <- c(min(fit$mu), max(fit$mu))

# compute zero inflation parameters omega

  if(is.matrix(Z)){ eta.omega <- Z %*% gamma }

  else { eta.omega <- Z * gamma }

  omega <- exp(eta.omega)/(1+exp(eta.omega))

  range.omega <- c(min(omega), max(omega))

# compute dispersion parameters phi

  if(is.matrix(W)){ eta.phi <- W %*% alpha }

  else { eta.phi <- W * alpha }

  phi <- 1+ exp(eta.phi)

  range.phi <- c(min(phi), max(phi))



# summarize results

  mle.data <- list(ZI.Parameter = omega,

    Coefficients.Mu = beta,

    Coefficients.Phi = alpha,

    Coefficients.Omega = gamma,

    Dispersion.Parameter = phi,

    Range.Mu = range.mu,

    Range.Phi = range.phi,

    Range.Omega = range.omega,

    Log.Likelihood = loglikelihood,

    Residuals = res,

    Pearson = chsq,

    AIC = AIC,

    Iterations = it,

    Time = time.needed,

    Message = message,

    Response = Y,

    Fitted.Values = fit$fit,

    Design.Mu = X,

    Design.Phi = W,

    Design.Omega = Z)

    if(summary) {

      summary.zigp(mle.data)

    }

    else {

      return(mle.data)

    }

  }

  else{

  return("X, W and Z have different numbers of rows!")

  }

}

