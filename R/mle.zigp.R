"mle.zigp" <-
function(Yin, Xin, Win=NULL, Zin=NULL, Offset = rep(1,length(Yin)), init = FALSE)

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

  if(is.null(W)==FALSE){
  
    if(is.matrix(W)) { k.alpha <<- dim(W)[2] }

    else { k.alpha <<- 1 }

  }

  if(is.null(Z)==FALSE){

    if(is.matrix(Z)) { k.gamma <<- dim(Z)[2] }

    else { k.gamma <<- 1 }
    
  }



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

    if(is.null(W)==FALSE){ alpha <- rep(-2.5,k.alpha) }
    
    else { alpha <- NULL }

    if(is.null(Z)==FALSE){ gamma <- rep(-2.5,k.gamma) }

    else { gamma <- NULL }

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

  coef <- delta

  beta <- coef[1 : k.beta]

  if(is.null(W)==FALSE){ alpha <- coef[(k.beta + 1) : (k.beta + k.alpha)] }
  
  else { alpha <- NULL }

  if(is.null(Z)==FALSE){ gamma <- coef[(k.beta + k.alpha + 1) : (k.beta + k.alpha + k.gamma)] }

  else { gamma <- NULL }



  fit <- fit.zigp(delta)

  res <- Y - fit$fit

  chsq <- sum((Y - fit$fit)^2 / ((1-fit$omega)*fit$mu*(fit$phi^2+fit$mu*fit$omega)))



  range.mu <- c(min(fit$mu), max(fit$mu))

# compute zero inflation parameters omega

  if(is.null(Z)==FALSE){

    if(is.matrix(Z)){ eta.omega <- Z %*% gamma }

    else { eta.omega <- Z * gamma }

    omega <- exp(eta.omega)/(1+exp(eta.omega))

  }
  
  else { omega <- rep(0,n) }

  range.omega <- c(min(omega), max(omega))

# compute dispersion parameters phi

  if(is.null(W)==FALSE){

    if(is.matrix(W)){ eta.phi <- W %*% alpha }

    else { eta.phi <- W * alpha }

    phi <- 1+ exp(eta.phi)

  }

  else { phi <- rep(1,n) }

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

    return(mle.data)

 }

